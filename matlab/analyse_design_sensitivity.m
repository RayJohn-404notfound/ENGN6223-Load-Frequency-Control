function summary = analyse_design_sensitivity()
%ANALYSE_DESIGN_SENSITIVITY Predeclared grid; original gains are not replaced.
% Run without arguments. Writes only results/design_sensitivity/.
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'results','design_sensitivity');
assert(isfile(fullfile(out,'protocol.json')),'Save the predeclared protocol before running.');
protocol=jsondecode(fileread(fullfile(out,'protocol.json')));assert(protocol.savedBeforeExecution);
src=jsondecode(fileread(fullfile(root,'results','research_summary.json')));p=src.p;spec=src.specifications;
old=readtable(fullfile(root,'results','parameter_sensitivity.csv'),'TextType','string');
base=old(old.controller=="PI",:);assert(height(base)==43);tau=(0:p.dt:p.tEnd-p.stepTime)';
[~,~,~,~,P]=plant(p);designCells=cell(0,1);id=0;
for wc=protocol.grid.wc_rad_s(:)'
 for pm=protocol.grid.targetPM_deg(:)'
  z=exp(1i*deg2rad(pm-180))/evalfr(P,1i*wc);
  id=id+1;designCells{id}=design('PI',wc,pm,0,real(z),-wc*imag(z),0,0);
  for N=protocol.grid.PIDF_N(:)'
   alpha=1/N;kp=(real(z)-alpha*imag(z))/(1+alpha);
   kd=(real(z)+imag(z))*(1+alpha^2)/(wc*(1+alpha));
   id=id+1;designCells{id}=design('PIDF',wc,pm,N,kp,wc*kp,kd,1/(N*wc));
  end
 end
end
designs=vertcat(designCells{:});assert(numel(designs)==36);rows=cell(0,1);crossRows=cell(0,1);checks=cell(0,1);
maxCrossResidual=0;maxModalLsimError=0;maxEigenCondition=0;maxEigenResidual=0;worstConditionCheck=struct;
for di=1:numel(designs)
 c=designs(di);
 if ~c.PositiveGains,continue;end
 fprintf('Design sensitivity %d/36: %s\n',di,c.DesignID);
 for k=1:height(base)
  pv=p;names={'Tg','Tt','Tp','Ksys','R'};
  for j=1:5,pv.(names{j})=p.(names{j})*base.([names{j} '_factor'])(k);end
  [A,B,E,C,Pv]=plant(pv);[Ac,Ec,Cy,Cu]=closed(A,B,E,C,c);
  lam=eig(Ac);stable=all(real(lam)<0);[num,den]=loopPolynomials(pv,c);
  fm=exactMargins(num,den);maxCrossResidual=max(maxCrossResidual,fm.maxResidual);
  r=blankRow(c,k,base);r.Stable=stable;r.MaxRealPole_per_s=max(real(lam));
  r.UnityCrossingCount=numel(fm.wpm);r.PM_deg=fm.pm;r.GM_dB=fm.gmDb;
  Ci=controller(c);L=Pv*Ci;
  if stable
   [V,D]=eig(Ac);lambda=diag(D);xss=-Ac\(Ec*p.loadStep);a=V\xss;
   r.EigenvectorCondition=cond(V);r.EigenResidual=norm(Ac*V-V*D,inf);
   maxEigenResidual=max(maxEigenResidual,r.EigenResidual);
   assert(r.EigenvectorCondition<1e8 && r.EigenResidual<1e-8);
   modes=exp(tau*lambda.');coef=[(Cy*V).'.*a,(Cu*V).'.*a];
   Y=real(repmat(([Cy;Cu]*xss).',numel(tau),1)-modes*coef);
   r.MaxAbs_Hz=max(abs(Y(:,1)));r.Nadir_Hz=min(Y(:,1));r.PeakCommand_pu=max(abs(Y(:,2)));
   last=find(abs(Y(:,1))>spec.recoveryBandHz,1,'last');
   if isempty(last),r.Recovery_s=0;elseif last<numel(tau),r.Recovery_s=tau(last+1);else,r.Recovery_s=NaN;end
   r.HorizonLimited=~isfinite(r.Recovery_s);r.IAE_Hz_s=trapz(tau,abs(Y(:,1)));
   r.TailBoundAtEndpoint_Hz=abs(Cy*xss)+sum(abs(coef(:,1)).*exp(real(lambda)*tau(end)));
   r.SteadyFrequency_Hz=Cy*xss;r.SteadyCommand_pu=Cu*xss;
   S=feedback(1,L);[r.Ms,r.MsFrequency_rad_s]=getPeakGain(S,1e-7);
   r.CommandNoiseGain100=abs(evalfr(Ci*S,100i));
   r.PassRecovery=r.Recovery_s<=spec.recoveryTime_s;
   r.PassDeviation=r.MaxAbs_Hz<=spec.maxAbsoluteDeviationHz;
   r.PassCommand=r.PeakCommand_pu<=spec.maxSecondaryCommandPu;
   r.PassPM=r.PM_deg>=spec.minPhaseMarginDeg;r.PassGM=r.GM_dB>=spec.minGainMarginDb;
   r.AllTargetsPass=stable && abs(r.SteadyFrequency_Hz)<1e-9 && r.PassRecovery && r.PassDeviation && r.PassCommand && r.PassPM && r.PassGM;
   if (c.IsOriginal && any(k==[1 28])) || r.EigenvectorCondition>maxEigenCondition
    Yls=lsim(ss(Ac,Ec,[Cy;Cu],zeros(2,1)),p.loadStep*ones(size(tau)),tau);
    err=max(abs(Y-Yls),[],'all');maxModalLsimError=max(maxModalLsimError,err);assert(err<1e-10);
    if r.EigenvectorCondition>maxEigenCondition
     maxEigenCondition=r.EigenvectorCondition;worstConditionCheck=struct('DesignID',c.DesignID,'CaseID',k,'EigenvectorCondition',maxEigenCondition,'ModalLsimDifference',err);
    end
   end
  end
  rows{end+1,1}=r; %#ok<AGROW>
  if k==1
   target=exp(1i*deg2rad(c.TargetPM_deg-180));pointError=abs(polyval(num,1i*c.Wc_rad_s)/polyval(den,1i*c.Wc_rad_s)-target);
   designs(di).NominalStable=stable;designs(di).PointMatchResidual=pointError;
   designs(di).UnityCrossings=numel(fm.wpm);designs(di).ActualPM_deg=fm.pm;
   designs(di).ActualGM_dB=fm.gmDb;
   match=any(abs(fm.wpm-c.Wc_rad_s)<=1e-7*c.Wc_rad_s);
   designs(di).MatchedEligible=stable && pointError<1e-10 && match && abs(fm.pm-c.TargetPM_deg)<1e-6;
   if ~stable,designs(di).Status='Excluded: nominal closed loop unstable';
   elseif ~designs(di).MatchedEligible,designs(di).Status='Unmatched: requested PM is not verified governing margin';
   else,designs(di).Status='Eligible';end
   am=allmargin(L);
   crossRows{end+1,1}=struct('DesignID',c.DesignID,'ExactGainCrossings_rad_s',jsonencode(fm.wpm), ...
    'ExactPhaseMargins_deg',jsonencode(fm.pms),'ExactPhaseCrossings_rad_s',jsonencode(fm.wgm), ...
    'ExactGainMargins_dB',jsonencode(fm.gmDbs),'MATLABGainCrossings_rad_s',jsonencode(am.PMFrequency), ...
    'MATLABPhaseMargins_deg',jsonencode(am.PhaseMargin),'MATLABPhaseCrossings_rad_s',jsonencode(am.GMFrequency), ...
    'MATLABGainMargins_dB',jsonencode(20*log10(am.GainMargin)),'PointMatchResidual',pointError); %#ok<AGROW>
  end
  if c.IsOriginal
   legacy=old(old.caseID==k & old.controller==string(c.Type),:);assert(height(legacy)==1);
   er=[r.MaxAbs_Hz-legacy.maxAbsHz,r.Recovery_s-legacy.recovery_s,r.IAE_Hz_s-legacy.IAE_Hz_s,r.PeakCommand_pu-legacy.uPeakPu];
   assert(all(abs(er)<1e-9),'Original time-domain benchmark no longer matches.');
   legacyPass=logical(legacy.stable) && legacy.recovery_s<=spec.recoveryTime_s && legacy.maxAbsHz<=spec.maxAbsoluteDeviationHz ...
    && legacy.uPeakPu<=spec.maxSecondaryCommandPu && legacy.PM_deg>=spec.minPhaseMarginDeg && legacy.GM_dB>=spec.minGainMarginDb;
   checks{end+1,1}=struct('DesignID',c.DesignID,'CaseID',k,'MaxAbsDifference_Hz',er(1),'RecoveryDifference_s',er(2), ...
    'IAEDifference_Hz_s',er(3),'CommandDifference_pu',er(4),'ExactMinusLegacyPM_deg',r.PM_deg-legacy.PM_deg, ...
    'ExactMinusLegacyGM_dB',r.GM_dB-legacy.GM_dB,'ClassificationUnchanged',r.AllTargetsPass==legacyPass); %#ok<AGROW>
  end
 end
 writetable(struct2table(designs),fullfile(out,'design_register.csv'));
 writetable(struct2table(vertcat(rows{:})),fullfile(out,'all_plant_cases.csv'));
end
T=struct2table(vertcat(rows{:}));D=struct2table(designs);R=struct2table(vertcat(checks{:}));
assert(all(R.ClassificationUnchanged));
writetable(struct2table(vertcat(crossRows{:})),fullfile(out,'nominal_crossing_inventory.csv'));
writetable(R,fullfile(out,'reference_checks.csv'));nom=T(T.CaseID==1,:);writetable(nom,fullfile(out,'nominal_metrics.csv'));
agg=cell(0,1);
for i=1:height(D)
 rr=T(strcmp(T.DesignID,D.DesignID{i}),:);n=height(rr);
 a=struct('DesignID',D.DesignID{i},'Type',D.Type{i},'Wc_rad_s',D.Wc_rad_s(i),'TargetPM_deg',D.TargetPM_deg(i), ...
  'FilterN',D.FilterN(i),'Status',D.Status{i},'IsOriginal',D.IsOriginal(i),'EvaluatedCases',n,'StableCases',0, ...
  'AllTargetsPassCases',0,'RecoveryPassCases',0,'DeviationPassCases',0,'PMPassCases',0,'GMPassCases',0, ...
  'HorizonLimitedCases',0,'WorstRecovery_s',NaN,'WorstMaxAbs_Hz',NaN,'WorstPM_deg',NaN,'WorstGM_dB',NaN,'WorstMs',NaN);
 if n>0
  a.StableCases=sum(rr.Stable);a.AllTargetsPassCases=sum(rr.AllTargetsPass);a.RecoveryPassCases=sum(rr.PassRecovery);
  a.DeviationPassCases=sum(rr.PassDeviation);a.PMPassCases=sum(rr.PassPM);a.GMPassCases=sum(rr.PassGM);a.HorizonLimitedCases=sum(rr.HorizonLimited);
  a.WorstRecovery_s=max(rr.Recovery_s,[],'omitnan');if any(rr.HorizonLimited),a.WorstRecovery_s=NaN;end;a.WorstMaxAbs_Hz=max(rr.MaxAbs_Hz,[],'omitnan');
  a.WorstPM_deg=min(rr.PM_deg);a.WorstGM_dB=min(rr.GM_dB);a.WorstMs=max(rr.Ms,[],'omitnan');
 end
 agg{end+1,1}=a; %#ok<AGROW>
end
AT=struct2table(vertcat(agg{:}));writetable(AT,fullfile(out,'design_summary.csv'));
pairs=cell(0,1);
for i=1:height(D)
 if ~strcmp(D.Type{i},'PIDF'),continue;end
 pi=find(strcmp(D.Type,'PI') & D.Wc_rad_s==D.Wc_rad_s(i) & D.TargetPM_deg==D.TargetPM_deg(i));assert(isscalar(pi));
 pair=struct('PI_DesignID',D.DesignID{pi},'PIDF_DesignID',D.DesignID{i},'Wc_rad_s',D.Wc_rad_s(i), ...
  'TargetPM_deg',D.TargetPM_deg(i),'FilterN',D.FilterN(i),'Eligible',D.MatchedEligible(pi)&&D.MatchedEligible(i), ...
  'Status','Not eligible: see both design-register statuses','PI_NominalRecovery_s',NaN,'PIDF_NominalRecovery_s',NaN, ...
  'PI_NominalIAE_Hz_s',NaN,'PIDF_NominalIAE_Hz_s',NaN,'PI_Noise100',NaN,'PIDF_Noise100',NaN, ...
  'PI_NominalMs',NaN,'PIDF_NominalMs',NaN,'PI_AllPassCases',NaN,'PIDF_AllPassCases',NaN, ...
  'PIDF_FasterCases',NaN,'PI_FasterCases',NaN,'EqualRecoveryCases',NaN,'ComparableRecoveryCases',NaN, ...
  'PIDF_LowerIAECases',NaN,'PIDF_LowerDeviationCases',NaN,'PIDF_HigherNoiseCases',NaN);
 if pair.Eligible
  a=T(strcmp(T.DesignID,pair.PI_DesignID),:);b=T(strcmp(T.DesignID,pair.PIDF_DesignID),:);assert(isequal(a.CaseID,b.CaseID));
  an=a(a.CaseID==1,:);bn=b(b.CaseID==1,:);valid=a.Stable & b.Stable & isfinite(a.Recovery_s) & isfinite(b.Recovery_s);
  pair.Status='Matched actual nominal crossover and governing PM';pair.PI_NominalRecovery_s=an.Recovery_s;pair.PIDF_NominalRecovery_s=bn.Recovery_s;
  pair.PI_NominalIAE_Hz_s=an.IAE_Hz_s;pair.PIDF_NominalIAE_Hz_s=bn.IAE_Hz_s;pair.PI_Noise100=an.CommandNoiseGain100;pair.PIDF_Noise100=bn.CommandNoiseGain100;
  pair.PI_NominalMs=an.Ms;pair.PIDF_NominalMs=bn.Ms;pair.PI_AllPassCases=sum(a.AllTargetsPass);pair.PIDF_AllPassCases=sum(b.AllTargetsPass);
  pair.PIDF_FasterCases=sum(valid & b.Recovery_s<a.Recovery_s-1e-9);pair.PI_FasterCases=sum(valid & a.Recovery_s<b.Recovery_s-1e-9);
  pair.EqualRecoveryCases=sum(valid & abs(a.Recovery_s-b.Recovery_s)<=1e-9);pair.ComparableRecoveryCases=sum(valid);
  pair.PIDF_LowerIAECases=sum(a.Stable & b.Stable & b.IAE_Hz_s<a.IAE_Hz_s);pair.PIDF_LowerDeviationCases=sum(a.Stable & b.Stable & b.MaxAbs_Hz<a.MaxAbs_Hz);
  pair.PIDF_HigherNoiseCases=sum(a.Stable & b.Stable & b.CommandNoiseGain100>a.CommandNoiseGain100);
 end
 pairs{end+1,1}=pair; %#ok<AGROW>
end
PT=struct2table(vertcat(pairs{:}));writetable(PT,fullfile(out,'matched_comparisons.csv'));
summary=struct('completedAt',char(datetime('now')),'matlabVersion',version,'protocolSavedAt',protocol.savedAt, ...
 'plannedPI',9,'plannedPIDF',27,'positivePI',sum(strcmp(D.Type,'PI') & D.PositiveGains), ...
 'positivePIDF',sum(strcmp(D.Type,'PIDF') & D.PositiveGains),'plannedMatchedPairs',27,'eligibleMatchedPairs',sum(PT.Eligible), ...
 'evaluatedPlantCases',height(T),'unstablePlantCases',sum(~T.Stable),'horizonLimitedCases',sum(T.HorizonLimited), ...
 'maxCrossingResidual',maxCrossResidual,'maxModalLsimDifference',maxModalLsimError, ...
 'maxEigenvectorCondition',maxEigenCondition,'maxEigenResidual',maxEigenResidual,'worstConditionCheck',worstConditionCheck, ...
 'maxReferenceResponseMetricDifference',max(abs(R{:,3:6}),[],'all'),'referenceClassificationsUnchanged',all(R.ClassificationUnchanged), ...
 'designs',table2struct(D),'aggregate',table2struct(AT),'matchedComparisons',table2struct(PT));
summary.scope='Finite exploratory design grid, conditional on the specified tuning closure and benchmark. No original-controller replacement or optimisation claim.';
summary.noiseBoundary='Secondary measurement only; command gain at 100 rad/s is one-frequency magnitude, not integrated noise variance.';
summary.filterBoundary='Changing N re-solves Kp,Ki,Kd, so N comparisons concern complete matched redesigns rather than an isolated filter effect.';
summary.exclusionBoundary='All 36 planned designs and27 planned pairs are retained. Infeasible candidates have zero evaluated cases and explicit exclusions, not fabricated simulation failures.';
writeJson(fullfile(out,'design_sensitivity_summary.json'),summary);
plot_design_sensitivity();writeFindings(out,summary,PT,AT,T);
summary=summarise_design_sensitivity();
fprintf('Design sensitivity finished: %d positive designs, %d cases, %d matched pairs.\n',sum(D.PositiveGains),height(T),sum(PT.Eligible));
end

function c=design(type,wc,pm,N,kp,ki,kd,tfilt)
if strcmp(type,'PI'),key=sprintf('PI_w%02d_pm%d',round(10*wc),pm);else,key=sprintf('PIDF_w%02d_pm%d_N%d',round(10*wc),pm,N);end
positive=all(isfinite([kp ki kd tfilt])) && kp>0 && ki>0 && (strcmp(type,'PI') || (kd>0 && tfilt>0));
status='Pending evaluation';if ~positive,status='Excluded: nonpositive required gain';end
c=struct('DesignID',key,'Type',type,'Wc_rad_s',wc,'TargetPM_deg',pm,'FilterN',N,'Kp',kp,'Ki',ki,'Kd',kd,'Tf',tfilt, ...
 'PositiveGains',positive,'NominalStable',false,'PointMatchResidual',NaN,'UnityCrossings',0,'ActualPM_deg',NaN,'ActualGM_dB',NaN, ...
 'MatchedEligible',false,'Status',status,'IsOriginal',wc==2.4 && pm==65 && (strcmp(type,'PI') || N==10));
end
function r=blankRow(c,k,b)
r=struct('DesignID',c.DesignID,'Type',c.Type,'Wc_rad_s',c.Wc_rad_s,'TargetPM_deg',c.TargetPM_deg,'FilterN',c.FilterN,'CaseID',b.caseID(k), ...
 'Tg_factor',b.Tg_factor(k),'Tt_factor',b.Tt_factor(k),'Tp_factor',b.Tp_factor(k),'Ksys_factor',b.Ksys_factor(k),'R_factor',b.R_factor(k), ...
 'Stable',false,'MaxRealPole_per_s',NaN,'UnityCrossingCount',0,'PM_deg',NaN,'GM_dB',NaN, ...
 'EigenvectorCondition',NaN,'EigenResidual',NaN,'MaxAbs_Hz',NaN,'Nadir_Hz',NaN,'Recovery_s',NaN,'HorizonLimited',true,'IAE_Hz_s',NaN,'PeakCommand_pu',NaN, ...
 'TailBoundAtEndpoint_Hz',NaN,'SteadyFrequency_Hz',NaN,'SteadyCommand_pu',NaN,'Ms',NaN,'MsFrequency_rad_s',NaN,'CommandNoiseGain100',NaN, ...
 'PassRecovery',false,'PassDeviation',false,'PassCommand',false,'PassPM',false,'PassGM',false,'AllTargetsPass',false);
end
function [A,B,E,C,P]=plant(p)
A=[-1/p.Tg 0 -1/(p.R*p.Tg);1/p.Tt -1/p.Tt 0;0 p.Ksys/p.Tp -1/p.Tp];B=[1/p.Tg;0;0];E=[0;0;-p.Ksys/p.Tp];C=[0 0 1];
P=tf(p.Ksys,[p.Tg*p.Tt*p.Tp,p.Tg*p.Tt+p.Tg*p.Tp+p.Tt*p.Tp,p.Tg+p.Tt+p.Tp,1+p.Ksys/p.R]);
end
function [Ac,Ec,Cy,Cu]=closed(A,B,E,C,c)
if strcmp(c.Type,'PI'),Ac=[A-B*c.Kp*C B*c.Ki;-C 0];Ec=[E;0];Cy=[C 0];Cu=[-c.Kp*C c.Ki];
else,q=c.Kd/c.Tf;Ac=[A-B*(c.Kp+q)*C B*c.Ki -B*q;-C 0 0;-C/c.Tf 0 -1/c.Tf];Ec=[E;0;0];Cy=[C 0 0];Cu=[-(c.Kp+q)*C c.Ki -q];end
end
function C=controller(c)
s=tf('s');C=c.Kp+c.Ki/s;if strcmp(c.Type,'PIDF'),C=C+c.Kd*s/(1+c.Tf*s);end
end
function [n,d]=loopPolynomials(p,c)
pd=[p.Tg*p.Tt*p.Tp,p.Tg*p.Tt+p.Tg*p.Tp+p.Tt*p.Tp,p.Tg+p.Tt+p.Tp,1+p.Ksys/p.R];
if strcmp(c.Type,'PI'),n=p.Ksys*[c.Kp c.Ki];d=conv([1 0],pd);
else,n=p.Ksys*[c.Kp*c.Tf+c.Kd,c.Kp+c.Ki*c.Tf,c.Ki];d=conv([c.Tf 1 0],pd);end
end
function m=exactMargins(n,d)
[ne,no]=split(n);[de,do]=split(d);
unity=add(add(conv(ne,ne),[0 conv(no,no)]),-add(conv(de,de),[0 conv(do,do)]));
q=positiveRoots(unity);wg=sqrt(q);z=polyval(n,1i*wg)./polyval(d,1i*wg);
pms=mod(rad2deg(angle(z))+360,360)-180; % equivalent signed PM in [-180,180)
res=abs(abs(z)-1);assert(all(res<1e-7));
phase=add(conv(no,de),-conv(ne,do));qp=positiveRoots(phase);wp=sqrt(qp);zp=polyval(n,1i*wp)./polyval(d,1i*wp);
keep=real(zp)<0 & isfinite(zp);wp=wp(keep);zp=zp(keep);gmDbs=20*log10(-1./real(zp));
phaseRes=abs(imag(zp))./max(1,abs(zp));assert(all(phaseRes<1e-7));
if isempty(pms),pm=Inf;else,[~,i]=min(abs(pms));pm=pms(i);end
if isempty(gmDbs),gmDb=Inf;else,[~,i]=min(abs(gmDbs));gmDb=gmDbs(i);end
m=struct('wpm',wg,'pms',pms,'wgm',wp,'gmDbs',gmDbs,'pm',pm,'gmDb',gmDb,'maxResidual',max([0;res(:);phaseRes(:)]));
end
function [e,o]=split(a)
a=fliplr(a);e=a(1:2:end).*(-1).^(0:numel(a(1:2:end))-1);o=a(2:2:end).*(-1).^(0:numel(a(2:2:end))-1);if isempty(o),o=0;end
end
function c=add(a,b)
c=zeros(1,max(numel(a),numel(b)));c(1:numel(a))=a;c(1:numel(b))=c(1:numel(b))+b;
end
function q=positiveRoots(a)
while numel(a)>1 && abs(a(end))<1e-14*max(abs(a)),a(end)=[];end
r=roots(fliplr(a));q=sort(real(r(abs(imag(r))<1e-7*max(1,abs(r)) & real(r)>1e-12)));
end
function writeJson(path,s)
f=fopen(path,'w');assert(f>=0);clean=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(s,'PrettyPrint',true));
end
function writeFindings(out,s,P,A,T)
f=fopen(fullfile(out,'findings.txt'),'w');assert(f>=0);clean=onCleanup(@()fclose(f));
fprintf(f,'EXPLORATORY DESIGN-SENSITIVITY STUDY\n\n');
fprintf(f,'The protocol was saved before simulation. All 36 planned candidates and 27 planned matched comparisons are retained.\n');
fprintf(f,'Positive-gain candidates: %d PI and %d PIDF. Eligible matched pairs: %d. Evaluated controller/plant cases: %d.\n',s.positivePI,s.positivePIDF,s.eligibleMatchedPairs,s.evaluatedPlantCases);
fprintf(f,'Unstable evaluated cases: %d. Horizon-limited rows: %d.\n',s.unstablePlantCases,s.horizonLimitedCases);
fprintf(f,'All nine PIDF designs at wc=1.8 rad/s require negative derivative gain. PI at wc=3.0 rad/s and PM=75 degrees requires negative integral gain. They are excluded, not replaced.\n');
e=P(P.Eligible,:);fprintf(f,'Across %d eligible matched nominal pairs, PIDF has shorter recovery in %d, lower IAE in %d, and higher one-frequency command-noise gain in %d.\n',height(e),sum(e.PIDF_NominalRecovery_s<e.PI_NominalRecovery_s),sum(e.PIDF_NominalIAE_Hz_s<e.PI_NominalIAE_Hz_s),sum(e.PIDF_Noise100>e.PI_Noise100));
fprintf(f,'These are conditional paired comparisons, not independent trials or probabilities. No pooled 8-versus-18 performance rate is claimed.\n');
fprintf(f,'Among positive candidates, PI all-target counts span %d to %d; PIDF counts span %d to %d. These ranges describe design sensitivity, not a fair pooled class comparison.\n',min(A.AllTargetsPassCases(strcmp(A.Type,'PI') & A.EvaluatedCases>0)),max(A.AllTargetsPassCases(strcmp(A.Type,'PI') & A.EvaluatedCases>0)),min(A.AllTargetsPassCases(strcmp(A.Type,'PIDF') & A.EvaluatedCases>0)),max(A.AllTargetsPassCases(strcmp(A.Type,'PIDF') & A.EvaluatedCases>0)));
fprintf(f,'Maximum original benchmark time-domain metric difference: %.5g; maximum modal/lsim response difference: %.5g. Original all-target classifications are unchanged.\n',s.maxReferenceResponseMetricDifference,s.maxModalLsimDifference);
fprintf(f,'Maximum complex crossing residual: %.5g. Complete nominal crossing inventories are exported; every parameter-case margin accounts for all crossings.\n',s.maxCrossingResidual);
fprintf(f,'Positive-gain cases with endpoint modal tail bound above the recovery band: %d. Sampled recovery remains explicitly horizon limited where applicable.\n',sum(T.Stable & T.TailBoundAtEndpoint_Hz>.001));
fprintf(f,'N changes Tf and the solved gains together. Results concern complete redesigns, not isolated derivative-filter causation.\n');
fprintf(f,'Noise enters the secondary measurement only. The 100 rad/s magnitude is not integrated noise variance or a sensor experiment.\n');
fprintf(f,'The retained original design remains the report benchmark. No optimisation, physical validation, or general superiority claim follows from this finite grid.\n');
end
