function summary = analyse_recovery_exceptions()
%ANALYSE_RECOVERY_EXCEPTIONS Explain selected fixed-gain recovery rankings.
% Run with no arguments from any working folder. Uses existing parameter and
% gain CSVs without changing them. Requires MATLAB and Control System Toolbox.
% All time coordinates below are elapsed time after the load step at t=1 s.
projectDir=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(projectDir,'results','deepening');if ~isfolder(out),mkdir(out);end
source=readtable(fullfile(projectDir,'results','parameter_sensitivity.csv'),'TextType','string');
gains=readtable(fullfile(projectDir,'results','controller_gains.csv'),'TextType','string');
caseIDs=[1 15 30 31];bands=[.0005 .0009 .001 .0011 .002];
names={'Tg','Tt','Tp','Ksys','R'};base=[.08 .3 20 120 2.4];
dt=.001;fineDt=.0001;absoluteEnd=60;loadTime=1;horizon=absoluteEnd-loadTime;
t=(0:dt:horizon)';tfine=(0:fineDt:horizon)';loadStep=.01;
recRows={};poleRows={};extRows={};crossRows={};freqRows={};validationRows={};
cases=cell(numel(caseIDs),1);plotData=cell(numel(caseIDs),2);allRootDifferences=[];
for ic=1:numel(caseIDs)
    caseID=caseIDs(ic);rr=source(source.caseID==caseID & source.controller=="PI",:);
    factors=rr{1,{'Tg_factor','Tt_factor','Tp_factor','Ksys_factor','R_factor'}};
    p=cell2struct(num2cell(base.*factors),names,2);
    caseResult=struct('caseID',caseID,'factorOrder',{names},'factors',factors,'parameters',p);
    controllerResults=cell(2,1);trajectory=table(t,t+loadTime,'VariableNames',{'Elapsed_s','AbsoluteTime_s'});
    modalTable=trajectory;
    for ci=1:2
        c=table2struct(gains(ci,:));control=char(c.name);
        [Ac,Ec,Cy,Cu,stepNum,stepDen,P,Cc]=model(p,c,loadStep);
        [V,D]=eig(Ac);lambda=diag(D);xss=-(Ac\Ec)*loadStep;
        alpha=V\xss;res=-(Cy*V).*alpha.';fss=Cy*xss;
        assert(all(real(lambda)<0),'Selected model is unstable.');
        exponent=exp(lambda*t');X=xss-V*(alpha.*exponent);Y=real([Cy;Cu]*X)';
        yf=real(fss+res*exp(lambda*tfine'));yf=yf(:);
        sys=ss(Ac,Ec,[Cy;Cu],zeros(2,1));ylsim=lsim(sys,loadStep*ones(size(t)),t);
        lsimError=max(abs(Y-ylsim),[],'all');
        checkTimes=unique([0 .01 .1 .7 1 2 3 4 5 8 15 30 horizon]);expmError=0;
        for te=checkTimes
            direct=[Cy;Cu]*(xss-expm(Ac*te)*xss);
            modal=real([Cy;Cu]*(xss-V*(alpha.*exp(lambda*te))));
            expmError=max(expmError,max(abs(direct-modal)));
        end
        [rt,pt]=residue(stepNum,stepDen);residueError=0;poleError=0;
        for k=1:numel(lambda)
            [pd,j]=min(abs(pt-lambda(k)));poleError=max(poleError,pd);
            residueError=max(residueError,abs(res(k)-rt(j)));
        end
        assert(lsimError<1e-9 && expmError<1e-10 && residueError<1e-8 && poleError<1e-7);
        envelope=@(q)abs(fss)+sum(abs(res(:)).*exp(real(lambda)*q));
        tailBound=envelope(horizon);assert(tailBound<min(bands));
        fun=@(q)real(fss+res*exp(lambda*q));
        deriv=@(q)real((res.*lambda.')*exp(lambda*q));
        d2=@(q)real((res.*(lambda.').^2)*exp(lambda*q));
        initialValueError=abs(fun(0));initialSlopeError=abs(deriv(0)-Cy*Ec*loadStep);
        assert(initialValueError<1e-12 && initialSlopeError<1e-12);
        [groupNames,groupMembers,components]=groupModes(lambda,res,t);
        % Extrema are found through derivative sign changes, then root refined.
        % The time cutoff excludes numerical zeros far below any tested band.
        extH=min(horizon,25);te=(0:dt:extH)';dv=real((res.*lambda.')*exp(lambda*te'))';
        extrema=bracketRoots(te,dv,deriv);extrema=extrema(extrema>1e-7);
        extremaData=zeros(numel(extrema),3);
        for ie=1:numel(extrema)
            et=extrema(ie);ev=fun(et);isMax=d2(et)<0;
            extremaData(ie,:)=[et ev isMax];
            extRows{end+1}=struct('CaseID',caseID,'Controller',control,'Elapsed_s',et, ...
                'Frequency_Hz',ev,'IsLocalMaximum',isMax); %#ok<AGROW>
        end
        reboundIndex=find(extremaData(:,1)>3 & extremaData(:,2)>0 & extremaData(:,3)==1,1);
        reboundTime=extremaData(reboundIndex,1);
        rebound=struct('elapsedTime_s',reboundTime,'frequency_Hz',fun(reboundTime), ...
            'excessAbove001Band_Hz',fun(reboundTime)-.001, ...
            'modeContributions_Hz',modalContributions(lambda,res,groupMembers,reboundTime));
        for k=1:numel(lambda)
            poleRows{end+1}=struct('CaseID',caseID,'Controller',control,'Mode',k, ...
                'PoleReal_per_s',real(lambda(k)),'PoleImag_rad_s',imag(lambda(k)), ...
                'StepResidueReal_Hz',real(res(k)),'StepResidueImag_Hz',imag(res(k)), ...
                'StepResidueMagnitude_Hz',abs(res(k))); %#ok<AGROW>
        end
        recoveryResults=cell(numel(bands),1);
        for ib=1:numel(bands)
            band=bands(ib);coarse=sampleRecovery(t,Y(:,1),band);fine=sampleRecovery(tfine,yf,band);
            crossing=[];fineCrossing=[];
            for side=[-1 1]
                rootFun=@(q)fun(q)-side*band;
                roots=bracketRoots(t,Y(:,1)-side*band,rootFun);
                rootsFine=bracketRoots(tfine,yf-side*band,rootFun);
                assert(numel(roots)==numel(rootsFine),'Grid refinement changed crossing count.');
                if ~isempty(roots),allRootDifferences=[allRootDifferences;abs(roots-rootsFine)];end %#ok<AGROW>
                crossing=[crossing;roots];fineCrossing=[fineCrossing;rootsFine]; %#ok<AGROW>
                for ir=1:numel(roots)
                    rt0=roots(ir);entering=side*deriv(rt0)<0;
                    crossRows{end+1}=struct('CaseID',caseID,'Controller',control,'Band_Hz',band, ...
                        'Elapsed_s',rt0,'BoundarySign',side,'EnteringBand',entering); %#ok<AGROW>
                end
            end
            crossing=sort(crossing);fineCrossing=sort(fineCrossing);
            if isempty(crossing),rootTime=0;else,rootTime=crossing(end);end
            assert(abs(coarse-rootTime)<=dt+1e-9 && abs(fine-rootTime)<=fineDt+1e-9);
            % Stable poles make this absolute modal envelope non-increasing.
            if envelope(0)<=band,certTime=0;else,certTime=fzero(@(q)envelope(q)-band,[0 horizon]);end
            rootContribution=modalContributions(lambda,res,groupMembers,rootTime);
            recoveryResults{ib}=struct('band_Hz',band,'coarseSampleRecovery_s',coarse, ...
                'fineSampleRecovery_s',fine,'rootRefinedRecovery_s',rootTime, ...
                'lastBoundarySign',sign(fun(rootTime)),'crossingTimes_s',crossing', ...
                'modeContributionsAtRecovery_Hz',rootContribution, ...
                'absoluteModalEnvelopeBelowBandAfter_s',certTime);
            recRows{end+1}=struct('CaseID',caseID,'Controller',control,'Band_Hz',band, ...
                'Recovery_dt001_s',coarse,'Recovery_dt0001_s',fine,'RecoveryRoot_s',rootTime, ...
                'CoarseMinusFine_s',coarse-fine,'LastBoundarySign',sign(fun(rootTime)), ...
                'CrossingCount',numel(crossing),'EnvelopeCertifiedAfter_s',certTime); %#ok<AGROW>
        end
        L=P*Cc;w=logspace(-3,3,10000);sv=abs(squeeze(freqresp(feedback(1,L),w)));
        [~,im]=max(sv);lo=log(w(max(1,im-1)));hi=log(w(min(numel(w),im+1)));
        wMs=exp(fminbnd(@(z)-abs(1/(1+evalfr(L,1i*exp(z)))),lo,hi));
        LM=evalfr(L,1i*wMs);Ms=1/abs(1+LM);[gm,pm,~,wcp]=margin(L);
        Ldesign=evalfr(L,2.4i);
        freqRows{end+1}=struct('CaseID',caseID,'Controller',control,'PM_deg',pm, ...
            'GM_dB',20*log10(gm),'GainCrossover_rad_s',wcp,'RefinedMs',Ms, ...
            'MsFrequency_rad_s',wMs,'LoopRealAtMs',real(LM),'LoopImagAtMs',imag(LM), ...
            'MinimumLocatedDistanceToMinusOne',abs(1+LM), ...
            'LoopMagnitudeAt2p4',abs(Ldesign),'LoopPhaseAt2p4_deg',rad2deg(angle(Ldesign)), ...
            'SensitivityMagnitudeAt2p4',1/abs(1+Ldesign)); %#ok<AGROW>
        validationRows{end+1}=struct('CaseID',caseID,'Controller',control, ...
            'ModalVsLsimMaxOutputError',lsimError,'ModalVsExpmMaxOutputError',expmError, ...
            'TFVsStatePoleMaxError',poleError,'TFVsStateResidueMaxError_Hz',residueError, ...
            'FrequencyEnvelopeAt60Absolute_s_Hz',tailBound,'EigenvectorCondition',cond(V), ...
            'InitialValueError_Hz',initialValueError,'InitialSlopeError_Hz_per_s',initialSlopeError); %#ok<AGROW>
        for ig=1:numel(groupNames)
            modalTable.(sprintf('%s_Mode%d_Hz',control,ig))=components(:,ig);
        end
        trajectory.([control '_Frequency_Hz'])=Y(:,1);trajectory.([control '_Command_pu'])=Y(:,2);
        controllerResults{ci}=struct('controller',control,'steadyFrequency_Hz',fss, ...
            'dominantPoleReal_per_s',max(real(lambda)),'polesRealImag',[real(lambda) imag(lambda)], ...
            'stepResiduesRealImag_Hz',[real(res(:)) imag(res(:))],'modeLabels',{groupNames}, ...
            'modeMembers',{groupMembers},'extremaTimeFrequencyIsMaximum',extremaData, ...
            'firstPositiveReboundAfter3s',rebound, ...
            'recovery',vertcat(recoveryResults{:}),'refinedMs',Ms,'MsFrequency_rad_s',wMs, ...
            'minimumLocatedDistanceToMinusOne',abs(1+LM),'modalTailBound_Hz',tailBound);
        plotData{ic,ci}=struct('t',t,'y',Y(:,1),'components',components,'modeLabels',{groupNames});
        old=source(source.caseID==caseID & source.controller==string(control),:);
        bandOne=recoveryResults{find(bands==.001,1)};
        assert(abs(old.recovery_s-bandOne.coarseSampleRecovery_s)<1e-9,'Existing recovery not reproduced.');
        fprintf('Case %d %s: dt=.001 %.3f; root %.9f s; dominant %.8f; Ms %.8f\n', ...
            caseID,control,bandOne.coarseSampleRecovery_s,bandOne.rootRefinedRecovery_s,max(real(lambda)),Ms);
    end
    caseResult.controllers=vertcat(controllerResults{:});
    switch caseID
        case 1
            caseResult.interpretation=['Nominal PIDF has faster observable oscillatory and slow real modes. ' ...
                'PI has a late negative trough near 4.38 s outside the band, whereas the PIDF later extrema remain inside 0.001 Hz. ' ...
                'Equal loop values at one crossover do not fix the complete disturbance-response poles and residues.'];
        case 15
            caseResult.interpretation=['PIDF has a slower observable oscillatory pair and later negative-trough recovery. ' ...
                'Its late positive rebound is just below 0.001 Hz, so narrowing the band to 0.0009 Hz creates an extra exit. ' ...
                'PIDF remains slower across the five tested bands.'];
        case 30
            caseResult.interpretation=['PIDF has faster dominant decay but less cancellation of its late positive oscillatory rebound by a negative real mode. ' ...
                'Its rebound exceeds 0.001 Hz by about 0.000056833 Hz, creating an additional exit absent for PI. ' ...
                'The ranking reverses at 0.0005, 0.0011 and 0.002 Hz, showing that the nominal-band exception is threshold-sensitive.'];
        case 31
            caseResult.interpretation=['PIDF has a slower observable oscillatory pair and a larger late positive rebound. ' ...
                'Both controllers leave the 0.001 Hz band during that rebound, but PIDF returns later. ' ...
                'PIDF remains slower across the five tested bands.'];
    end
    cases{ic}=caseResult;
    writetable(trajectory,fullfile(out,sprintf('exception_trajectory_case%02d.csv',caseID)));
    writetable(modalTable,fullfile(out,sprintf('exception_modal_case%02d.csv',caseID)));
end
recoveryTable=struct2table(vertcat(recRows{:}));
writetable(recoveryTable,fullfile(out,'exception_band_recovery.csv'));
writetable(struct2table(vertcat(poleRows{:})),fullfile(out,'exception_poles_residues.csv'));
writetable(struct2table(vertcat(extRows{:})),fullfile(out,'exception_extrema.csv'));
writetable(struct2table(vertcat(crossRows{:})),fullfile(out,'exception_band_crossings.csv'));
writetable(struct2table(vertcat(freqRows{:})),fullfile(out,'exception_frequency_metrics.csv'));
writetable(struct2table(vertcat(validationRows{:})),fullfile(out,'exception_validation.csv'));
summary=struct('method','Stable augmented state-space modal step response; residue, lsim and expm cross-checks', ...
    'absoluteLoadStepTime_s',loadTime,'absoluteEndTime_s',absoluteEnd,'elapsedHorizon_s',horizon, ...
    'loadStep_pu',loadStep,'coarse_dt_s',dt,'fine_dt_s',fineDt,'bands_Hz',bands, ...
    'rootGridRefinementMaxDifference_s',max(allRootDifferences),'caseResults',vertcat(cases{:}), ...
    'validation',vertcat(validationRows{:}),'frequencyMetrics',vertcat(freqRows{:}), ...
    'limitations','Selected deterministic benchmark cases; complete designs have different PI gains. Modal and threshold mechanisms do not isolate derivative action or establish physical causation.');
fid=fopen(fullfile(out,'exception_summary.json'),'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));clear clean
makeFigures(out,caseIDs,bands,recoveryTable,plotData);
writeFindings(out,summary,recoveryTable);
fprintf('Exception analysis complete: %s\n',out);
end

function [Ac,Ec,Cy,Cu,num,den,P,Cc]=model(p,c,d)
A=[-1/p.Tg 0 -1/(p.R*p.Tg);1/p.Tt -1/p.Tt 0;0 p.Ksys/p.Tp -1/p.Tp];
B=[1/p.Tg;0;0];E=[0;0;-p.Ksys/p.Tp];C=[0 0 1];
Dp=[p.Tg*p.Tt*p.Tp,p.Tg*p.Tt+p.Tg*p.Tp+p.Tt*p.Tp,p.Tg+p.Tt+p.Tp,1+p.Ksys/p.R];
P=tf(p.Ksys,Dp);N=conv([p.Tg 1],[p.Tt 1]);
if c.Kd==0
    Ac=[A-B*c.Kp*C B*c.Ki;-C 0];Ec=[E;0];Cy=[C 0];Cu=[-c.Kp*C c.Ki];
    den=polyadd(conv([1 0],Dp),p.Ksys*[c.Kp c.Ki]);num=-d*p.Ksys*N;
    Cc=pid(c.Kp,c.Ki);
else
    kd=c.Kd/c.Tf;
    Ac=[A-B*(c.Kp+kd)*C B*c.Ki -B*kd;-C 0 0;-C/c.Tf 0 -1/c.Tf];
    Ec=[E;0;0];Cy=[C 0 0];Cu=[-(c.Kp+kd)*C c.Ki -kd];
    den=polyadd(conv([c.Tf 1 0],Dp),p.Ksys*[c.Kp*c.Tf+c.Kd c.Kp+c.Ki*c.Tf c.Ki]);
    num=-d*p.Ksys*conv(N,[c.Tf 1]);Cc=pid(c.Kp,c.Ki,c.Kd,c.Tf);
end
end

function c=polyadd(a,b)
n=max(numel(a),numel(b));c=[zeros(1,n-numel(a)) a]+[zeros(1,n-numel(b)) b];
end

function roots=bracketRoots(t,y,fun)
idx=find(y(1:end-1).*y(2:end)<0);roots=zeros(numel(idx),1);
for k=1:numel(idx),roots(k)=fzero(fun,[t(idx(k)) t(idx(k)+1)],optimset('TolX',1e-12));end
exact=find(y==0);roots=unique([roots;t(exact)]);roots=sort(roots);
end

function value=sampleRecovery(t,y,band)
idx=find(abs(y)>band,1,'last');
if isempty(idx),value=0;elseif idx==numel(t),value=NaN;else,value=t(idx+1);end
end

function [labels,members,components]=groupModes(lambda,res,t)
members={};labels={};used=false(size(lambda));
[~,order]=sort(real(lambda),'descend');
for k=order(:)'
    if used(k),continue,end
    if abs(imag(lambda(k)))<1e-9
        ids=k;labels{end+1}=sprintf('Real pole %.4f',real(lambda(k))); %#ok<AGROW>
    else
        [~,j]=min(abs(lambda-conj(lambda(k))));ids=[k j];
        labels{end+1}=sprintf('Pair %.4f +/- j%.4f',real(lambda(k)),abs(imag(lambda(k)))); %#ok<AGROW>
    end
    used(ids)=true;members{end+1}=ids; %#ok<AGROW>
end
components=zeros(numel(t),numel(members));
for k=1:numel(members),ids=members{k};components(:,k)=real(res(ids)*exp(lambda(ids)*t'))';end
end

function values=modalContributions(lambda,res,members,t)
values=zeros(1,numel(members));
for k=1:numel(members),ids=members{k};values(k)=real(res(ids)*exp(lambda(ids)*t));end
end

function makeFigures(out,ids,bands,rt,pd)
col=[0 .35 .65;.85 .30 .08];
f=figure('Visible','off','Color','w','Position',[100 100 1150 730]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for k=1:numel(ids)
    nexttile;hold on
    for ci=1:2,plot(pd{k,ci}.t,pd{k,ci}.y,'LineWidth',1.5,'Color',col(ci,:));end
    yline(.001,'k:');yline(-.001,'k:');xlim([1.5 6]);ylim([-.0035 .006]);grid on
    xlabel('Time after load increase (s)');ylabel('Frequency deviation (Hz)');
    title(sprintf('Case %d: recovery-band detail',ids(k)));legend('PI','PIDF','Location','northeast');
end
saveFigure(f,out,'exception_recovery_detail');
f=figure('Visible','off','Color','w','Position',[100 100 1150 730]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for k=1:numel(ids)
    nexttile;hold on
    for ci=1:2
        cn=["PI" "PIDF"];r=rt(rt.CaseID==ids(k) & string(rt.Controller)==cn(ci),:);
        plot(r.Band_Hz*1000,r.RecoveryRoot_s,'-o','LineWidth',1.5,'Color',col(ci,:));
    end
    xlim([min(bands) max(bands)]*1000);ylim([0 8]);grid on
    xlabel('Absolute recovery band (mHz)');ylabel('Root-refined recovery (s)');
    title(sprintf('Case %d: threshold dependence',ids(k)));legend('PI','PIDF','Location','northeast');
end
saveFigure(f,out,'exception_band_sensitivity');
f=figure('Visible','off','Color','w','Position',[100 100 1250 1100]);
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
for k=1:numel(ids)
    for ci=1:2
        nexttile;data=pd{k,ci};plot(data.t,data.components,'LineWidth',1);hold on
        plot(data.t,data.y,'k','LineWidth',1.6);yline(0,'k:');grid on;xlim([0 8]);ylim([-.035 .035]);
        xlabel('Time after load increase (s)');ylabel('Frequency contribution (Hz)');
        cn={'PI','PIDF'};title(sprintf('Case %d, %s',ids(k),cn{ci}));
        legend([data.modeLabels {'Sum'}],'Location','northeast','FontSize',7);
    end
end
saveFigure(f,out,'exception_modal_decomposition');
end

function saveFigure(f,out,name)
style_lfc_figure(f);
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector');
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
savefig(f,fullfile(out,[name '.fig']));close(f)
end

function writeFindings(out,summary,rt)
fid=fopen(fullfile(out,'exception_findings.txt'),'w');assert(fid>=0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'RECOVERY EXCEPTION ANALYSIS\n\n');
fprintf(fid,'Method and evidence boundary\n');
fprintf(fid,'The gains and selected parameter cases are read from the existing project CSVs. No gain is retuned.\n');
fprintf(fid,'A 0.01 pu load step starts at absolute time 1 s. The selected simulations end at 60 s.\n');
fprintf(fid,'All recovery and extremum times below are measured from the load increase.\n');
fprintf(fid,'The response is f(t)=f_ss+sum r_k exp(lambda_k t), with r_k=-(Cy V)_k (V\\x_ss)_k.\n');
fprintf(fid,'The residues multiply exponentials and have units Hz. Complex conjugate modes are combined into real contributions.\n');
fprintf(fid,'Transfer-function residues, direct matrix exponentials and lsim independently check this modal representation.\n');
fprintf(fid,'Both 0.001 s and 0.0001 s grids bracket the same band crossings, followed by scalar root refinement.\n');
fprintf(fid,'The maximum root difference between these grids is %.3g s.\n',summary.rootGridRefinementMaxDifference_s);
fprintf(fid,'Matching crossing counts supports numerical convergence; it is not a formal proof excluding every possible unsampled tangency.\n');
fprintf(fid,'An absolute modal envelope bounds the remaining tail below every tested band after the simulated horizon.\n\n');
fprintf(fid,'This mathematical envelope is evaluated in floating point; it is not an interval-arithmetic certificate.\n\n');
fprintf(fid,'Root-refined recovery times (PI / PIDF, seconds)\n');
for k=1:numel(summary.caseResults)
    cr=summary.caseResults(k);fprintf(fid,'Case %d; parameter factors: %s\n',cr.caseID,mat2str(cr.factors,3));
    fprintf(fid,'  Interpretation: %s\n',cr.interpretation);
    for b=summary.bands_Hz
        r=rt(rt.CaseID==cr.caseID & abs(rt.Band_Hz-b)<1e-12,:);
        fprintf(fid,'  Band %.4g Hz: %.9f / %.9f\n',b,r.RecoveryRoot_s(1),r.RecoveryRoot_s(2));
    end
    for ci=1:2
        c=cr.controllers(ci);fprintf(fid,'  %s dominant pole real part: %.9f /s; Ms %.7f at %.7f rad/s.\n', ...
            c.controller,c.dominantPoleReal_per_s,c.refinedMs,c.MsFrequency_rad_s);
        ex=c.extremaTimeFrequencyIsMaximum;ix=find(ex(:,1)>1 & ex(:,1)<8 & abs(ex(:,2))>.00025);
        for j=ix'
            fprintf(fid,'    Extremum: t=%.9f s, f=%.9f Hz, maximum=%d.\n',ex(j,1),ex(j,2),ex(j,3));
        end
        one=c.recovery(abs([c.recovery.band_Hz]-.001)<1e-12);
        fprintf(fid,'    Modal contributions at final 0.001 Hz crossing: %s Hz.\n',mat2str(one.modeContributionsAtRecovery_Hz,9));
        fprintf(fid,'    Mode order: %s\n',strjoin(c.modeLabels,'; '));
        rb=c.firstPositiveReboundAfter3s;
        fprintf(fid,'    Late positive rebound: %.9f Hz at %.9f s; contributions %s Hz.\n', ...
            rb.frequency_Hz,rb.elapsedTime_s,mat2str(rb.modeContributions_Hz,9));
    end
end
fprintf(fid,'\nInterpretation boundary\n');
fprintf(fid,'Dominant poles describe asymptotic decay, while residues, phase and the recovery threshold determine the final band crossing.\n');
fprintf(fid,'Each response compares a complete design. PI and PIDF also have different proportional and integral gains.\n');
fprintf(fid,'The numerical findings do not isolate the causal contribution of derivative action or establish behaviour of a physical grid.\n');
fprintf(fid,'Equal loop magnitude and phase at one frequency do not determine the full return ratio or disturbance response.\n');
f=summary.frequencyMetrics;
fprintf(fid,'At nominal 2.4 rad/s, both loops have magnitude one, phase -115 degrees and sensitivity magnitude %.9f.\n',f(1).SensitivityMagnitudeAt2p4);
fprintf(fid,'The located nominal sensitivity peaks occur elsewhere: PI %.9f at %.9f rad/s; PIDF %.9f at %.9f rad/s.\n', ...
    f(1).RefinedMs,f(1).MsFrequency_rad_s,f(2).RefinedMs,f(2).MsFrequency_rad_s);
fprintf(fid,'The corresponding distances from -1 are %.9f and %.9f, rather than a common distance fixed by phase margin.\n', ...
    f(1).MinimumLocatedDistanceToMinusOne,f(2).MinimumLocatedDistanceToMinusOne);
fprintf(fid,'The refined sensitivity peaks are local refinements around the largest sampled peaks, not certified global maxima.\n');
fprintf(fid,'The exported extrema table is limited to the first 25 seconds; recovery searches cover the full 59-second post-step horizon.\n');
fprintf(fid,'The band figure connects five tested widths. These segments guide the eye and do not assert continuous interpolation between samples.\n');
end
