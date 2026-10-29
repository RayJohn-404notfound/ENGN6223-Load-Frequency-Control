%% ENGN6223: single-area load frequency control study
% Run this script from any folder. All outputs are saved beside this project.
% Requires MATLAB and Control System Toolbox. Simulink cross-check is optional.
% The original baseline folder is not modified.
projectDir=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(projectDir,'matlab'));
out=fullfile(projectDir,'results'); if ~isfolder(out),mkdir(out);end
if ~isfolder(fullfile(projectDir,'model')),mkdir(fullfile(projectDir,'model'));end
diary(fullfile(out,'matlab_validation_log.txt'));
fprintf('LFC research | %s | MATLAB %s\n',char(datetime('now')),version);
assert(~isempty(which('ss')),'Control System Toolbox is required for this study.');
p=struct('Tg',.08,'Tt',.3,'Tp',20,'Ksys',120,'R',2.4, ...
    'loadStep',.01,'stepTime',1,'tEnd',30,'dt',.001);
[A,B,E,C,P,Q]=makePlant(p); s=tf('s');
rankCo=rank(ctrb(A,B)); rankOb=rank(obsv(A,C));
assert(rankCo==3 && rankOb==3); assert(all(real(eig(A))<0));
% P is the secondary-command-to-frequency plant with primary droop closed.
% Q is the load-disturbance-to-frequency transfer function, with its minus sign.
den=[p.Tg*p.Tt*p.Tp,p.Tg*p.Tt+p.Tg*p.Tp+p.Tt*p.Tp,p.Tg+p.Tt+p.Tp,1+p.Ksys/p.R];
Pformula=tf(p.Ksys,den); Qformula=-tf(p.Ksys*conv([p.Tg 1],[p.Tt 1]),den);
w=logspace(-3,3,2000)';
assert(max(abs(squeeze(freqresp(P-Pformula,w))))<1e-10);
assert(max(abs(squeeze(freqresp(Q-Qformula,w))))<1e-10);

%% Transparent loop shaping: same desired crossover and phase margin
wc=2.4; targetPM=65; Tf=1/(10*wc);
Cstar=exp(1i*deg2rad(-180+targetPM))/evalfr(P,1i*wc);
controllers(1)=struct('name','PI','Kp',real(Cstar),'Ki',-wc*imag(Cstar),'Kd',0,'Tf',0);
Dstar=1i*wc/(1+Tf*1i*wc);
gains=[1 real(Dstar);-1 imag(Dstar)]\[real(Cstar);imag(Cstar)];
controllers(2)=struct('name','PIDF','Kp',gains(1),'Ki',wc*gains(1),'Kd',gains(2),'Tf',Tf);
assert(all([controllers.Kp controllers.Ki]>0) && controllers(2).Kd>0);
% These are project specifications, not grid-code limits.
spec=struct('recoveryBandHz',.001,'recoveryTime_s',8,'maxAbsoluteDeviationHz',.025, ...
    'maxSecondaryCommandPu',.02,'minPhaseMarginDeg',45,'minGainMarginDb',6);

%% Nominal time responses; disturbance step is handled exactly
nominal=struct([]); loops=cell(1,2); cls=cell(1,3); frequencyData=[]; validation=struct;
for i=1:3
    if i==1,c=[];name='Droop';else,c=controllers(i-1);name=c.name;end
    [sys,Ac,Ec,Cy,Cu]=closedModel(p,c);cls{i}=sys;
    [t,Y,X]=response(sys,p);frequencyData(:,i)=Y(:,1);
    met=measures(t,Y,p,spec.recoveryBandHz);
    met.name=name;met.maxRealPole=max(real(eig(Ac)));met.poles=[real(eig(Ac)),imag(eig(Ac))];
    met.steadyFrequencyHz=-(Cy/Ac)*Ec*p.loadStep;
    met.steadySecondaryPu=-(Cu/Ac)*Ec*p.loadStep;
    met.PM_deg=NaN;met.GM_dB=NaN;met.wcg_rad_s=NaN;met.wcp_rad_s=NaN;
    met.Ms=NaN;met.Mt=NaN;met.noiseCommandGain100=NaN;
    if i>1
        Ci=controllerTF(c);L=minreal(P*Ci,1e-9);loops{i-1}=L;
        assert(abs(evalfr(L,1i*wc)-exp(1i*deg2rad(-180+targetPM)))<1e-10);
        [gm,pm,wgm,wpm]=margin(L);am=allmargin(L);
        met.PM_deg=pm;met.GM_dB=20*log10(gm);met.wcg_rad_s=wgm;met.wcp_rad_s=wpm;
        met.allMargins=am;
        S=feedback(1,L);T=feedback(L,1);CS=minreal(Ci*S,1e-9);
        met.Ms=max(abs(squeeze(freqresp(S,w))));met.Mt=max(abs(squeeze(freqresp(T,w))));
        met.noiseCommandGain100=abs(evalfr(CS,100i));
        ny=nyquistCheck(L);met.nyquist=ny;
        assert(ny.closedUnstablePoles==0 && ny.windingCCW==0 && ny.windingRefinedCCW==0);
        assert(abs(met.steadyFrequencyHz)<1e-10 && abs(met.steadySecondaryPu-p.loadStep)<1e-10);
        met.nominalTargetsPass=met.recovery_s<=spec.recoveryTime_s && met.maxAbsHz<=spec.maxAbsoluteDeviationHz ...
            && met.uPeakPu<=spec.maxSecondaryCommandPu && pm>=spec.minPhaseMarginDeg && met.GM_dB>=spec.minGainMarginDb;
        assert(met.nominalTargetsPass);
        % Transfer-function interconnection and explicit state equations must agree.
        outTF=minreal([Q/(1+P*Ci);-Ci*Q/(1+P*Ci)],1e-8);
        [~,Ytf]=response(ss(outTF),p);
        validation.([name '_transferStateMaxError'])=max(abs(Ytf(:)-Y(:)));
        assert(validation.([name '_transferStateMaxError'])<1e-7);
    else
        met.allMargins=struct;met.nyquist=struct;met.nominalTargetsPass=false;
        assert(abs(met.steadyFrequencyHz+.01/(1/120+1/2.4))<1e-10);
    end
    % Independently integrate the physical state equations after the step.
    ta=(p.stepTime:p.dt:p.tEnd)';
    [~,xo]=ode45(@(~,x)Ac*x+Ec*p.loadStep,ta,zeros(size(Ac,1),1), ...
        odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.01));
    idx=t>=p.stepTime;
    validation.([name '_odeLsimMaxStateError'])=max(abs(xo-X(idx,:)),[],'all');
    assert(validation.([name '_odeLsimMaxStateError'])<1e-7);
    yi=[Cy;Cu]*xo';
    nominal(i).name=name;nominal(i).metrics=met;nominal(i).Y=Y;nominal(i).X=X;
    fprintf('%s: nadir %.8f Hz, recovery %.3f s, upeak %.6f pu, PM %.2f deg, GM %.2f dB\n', ...
        name,met.nadirHz,met.recovery_s,met.uPeakPu,met.PM_deg,met.GM_dB);
end
nomTable=metricsTable(nominal);writetable(nomTable,fullfile(out,'nominal_metrics.csv'));
writetable(table(t,.01*double(t>=1),nominal(1).Y(:,1),nominal(2).Y(:,1),nominal(3).Y(:,1), ...
    nominal(2).Y(:,2),nominal(3).Y(:,2),'VariableNames', ...
    {'Time_s','Load_pu','Droop_df_Hz','PI_df_Hz','PIDF_df_Hz','PI_u_pu','PIDF_u_pu'}),fullfile(out,'nominal_response.csv'));
writetable(struct2table(controllers),fullfile(out,'controller_gains.csv'));
% PI quartic Routh first column, including its integrator.
cpi=controllers(1);coeff=[den(1:3),den(4)+p.Ksys*cpi.Kp,p.Ksys*cpi.Ki];
b1=(coeff(2)*coeff(3)-coeff(1)*coeff(4))/coeff(2);
c1=(b1*coeff(4)-coeff(2)*coeff(5))/b1;
routhFirstColumn=[coeff(1);coeff(2);b1;c1;coeff(5)];assert(all(routhFirstColumn>0));

%% Fixed-gain sensitivity: 32 corners, 10 axial cases and nominal
names={'Tg','Tt','Tp','Ksys','R'};corners=.8+.4*(dec2bin(0:31,5)-'0');
axial=ones(10,5);for j=1:5,axial(2*j-1,j)=.8;axial(2*j,j)=1.2;end
factors=[ones(1,5);corners;axial];robustRows=cell(size(factors,1)*2,1);envelopeMin=inf(numel(t),2);envelopeMax=-inf(numel(t),2);
ri=0;
for k=1:size(factors,1)
    pv=p;for j=1:5,pv.(names{j})=p.(names{j})*factors(k,j);end
    [~,~,~,~,Pv]=makePlant(pv);
    for ci=1:2
        ri=ri+1;[sysv,Av]=closedModel(pv,controllers(ci));stable=all(real(eig(Av))<0);
        [~,Yv]=response(sysv,pv);mv=measures(t,Yv,pv,spec.recoveryBandHz);
        [gm,pm]=margin(Pv*controllerTF(controllers(ci)));
        robustRows{ri}=struct('caseID',k,'controller',controllers(ci).name,'Tg_factor',factors(k,1), ...
            'Tt_factor',factors(k,2),'Tp_factor',factors(k,3),'Ksys_factor',factors(k,4),'R_factor',factors(k,5), ...
            'stable',stable,'maxRealPole',max(real(eig(Av))),'nadirHz',mv.nadirHz,'maxAbsHz',mv.maxAbsHz, ...
            'recovery_s',mv.recovery_s,'IAE_Hz_s',mv.IAE_Hz_s,'uPeakPu',mv.uPeakPu,'PM_deg',pm,'GM_dB',20*log10(gm));
        envelopeMin(:,ci)=min(envelopeMin(:,ci),Yv(:,1));envelopeMax(:,ci)=max(envelopeMax(:,ci),Yv(:,1));
    end
end
robustTable=struct2table(vertcat(robustRows{:}));writetable(robustTable,fullfile(out,'parameter_sensitivity.csv'));
robustSummaryRows=cell(2,1);
for ci=1:2
    rr=robustTable(strcmp(robustTable.controller,controllers(ci).name),:);
    robustSummaryRows{ci}=struct('name',controllers(ci).name,'testedCases',height(rr),'stableCases',sum(rr.stable), ...
        'worstNadirHz',min(rr.nadirHz),'worstRecovery_s',max(rr.recovery_s),'maxCommandPu',max(rr.uPeakPu), ...
        'minPM_deg',min(rr.PM_deg),'minGM_dB',min(rr.GM_dB), ...
        'nominalTimeTargetCases',sum(rr.recovery_s<=spec.recoveryTime_s), ...
        'allNominalTargetsCases',sum(logical(rr.stable) & rr.recovery_s<=spec.recoveryTime_s ...
            & rr.maxAbsHz<=spec.maxAbsoluteDeviationHz & rr.uPeakPu<=spec.maxSecondaryCommandPu ...
            & rr.PM_deg>=spec.minPhaseMarginDeg & rr.GM_dB>=spec.minGainMarginDb));
end
robustSummary=vertcat(robustSummaryRows{:});
fprintf('Parameter sensitivity: %d fixed-gain cases per controller.\n',size(factors,1));disp(struct2table(robustSummary));

%% Nonlinear actuator saturation, with and without anti-windup
% A project stress test, not a validated model of a specific actuator.
satSettings=struct('loadPu',.03,'on_s',5,'off_s',20,'end_s',50,'limitPu',.015,'Taw_s',.5,'dt',.005);
saturationRows=cell(4,1);n=0;
for ci=1:2
    for aw=[false,true]
        n=n+1;[ts,Ys,Zs]=saturationResponse(p,controllers(ci),satSettings,aw);
        after=ts>=satSettings.off_s;last=find(abs(Ys(after,1))>spec.recoveryBandHz,1,'last');ta=ts(after);
        rec=NaN;if isempty(last),rec=0;elseif last<numel(ta),rec=ta(last+1)-satSettings.off_s;end
        saturationRows{n}=struct('controller',controllers(ci).name,'antiWindup',aw,'t',ts,'Y',Ys, ...
            'integralOutput',Zs(:,4),'recoveryAfterRemoval_s',rec,'peakFrequencyHz',max(Ys(:,1)), ...
            'maxAbsHz',max(abs(Ys(:,1))),'maxCommandPu',max(abs(Ys(:,2))));
        assert(max(abs(Ys(:,2)))<=satSettings.limitPu+1e-12);
    end
end
saturation=vertcat(saturationRows{:});
satTable=table(string({saturation.controller})',[saturation.antiWindup]',[saturation.recoveryAfterRemoval_s]', ...
    [saturation.peakFrequencyHz]',[saturation.maxAbsHz]',[saturation.maxCommandPu]', ...
    'VariableNames',{'Controller','AntiWindup','RecoveryAfterRemoval_s','PositivePeak_Hz','MaxAbsoluteDeviation_Hz','MaxCommand_pu'});
writetable(satTable,fullfile(out,'saturation_metrics.csv'));
for n=1:4
    writetable(table(saturation(n).t,saturation(n).Y(:,1),saturation(n).Y(:,2),saturation(n).integralOutput, ...
        'VariableNames',{'Time_s','FrequencyDeviation_Hz','Command_pu','IntegralOutput_pu'}), ...
        fullfile(out,sprintf('saturation_%s_AW%d.csv',saturation(n).controller,saturation(n).antiWindup)));
end

%% Optional course extension: state feedback, integral action and observer
Aa=[A zeros(3,1);-C 0];Ba=[B;0];Ea=[E;0];rankAug=rank(ctrb(Aa,Ba));assert(rankAug==4);
assignedPoles=[-1.5 -2 -3 -8];Kplace=place(Aa,Ba,assignedPoles);
AclPlace=Aa-Ba*Kplace;
[~,Yplace]=response(ss(AclPlace,Ea,[[C 0];-Kplace],zeros(2,1)),p);
% Bryson-style scaling: x1,x2 and integral error have stated illustrative scales.
Qcost=diag(1./[.02 .02 .025 .025].^2);Rcost=1/.02^2;
[Klqr,~,lqrPoles]=lqr(Aa,Ba,Qcost,Rcost);
[~,Ylqr]=response(ss(Aa-Ba*Klqr,Ea,[[C 0];-Klqr],zeros(2,1)),p);
% Observer demonstration is disturbance-free, with an initial estimation error.
% It is not required by the implemented PI/PID controllers, which measure frequency.
observerPoles=[-6 -7 -9];Lobs=place(A',C',observerPoles)';
tobs=(0:.005:5)';x0=[0;0;.02];xhat0=zeros(3,1);
Aobs=[A zeros(3);Lobs*C A-Lobs*C];Zobs=initial(ss(Aobs,zeros(6,1),eye(6),zeros(6,1)),[x0;xhat0],tobs);
obsError=Zobs(:,1:3)-Zobs(:,4:6);
assert(norm(obsError(end,:))<1e-8);
separationPoles=eig([AclPlace,[B*Kplace(1:3);zeros(1,3)];zeros(3,4),A-Lobs*C]);
assert(max(abs(sort(real(separationPoles))-sort([assignedPoles observerPoles]')))<1e-7);
extension=struct('rankControllability',rankCo,'rankObservability',rankOb,'rankIntegralControllability',rankAug, ...
    'stateFeedbackK',Kplace,'assignedPoles',assignedPoles,'LQR_K',Klqr,'LQR_Q',Qcost,'LQR_R',Rcost, ...
    'LQR_poles',[real(lqrPoles),imag(lqrPoles)],'observerGain',Lobs,'observerPoles',observerPoles, ...
    'observerFinalErrorNorm',norm(obsError(end,:)),'observerDisturbanceKnownZero',true, ...
    'separationPoles',[real(separationPoles),imag(separationPoles)]);

%% Independent Simulink diagram cross-check
validation.Simulink='NOT_RUN';
if ~isempty(which('sim')) && license('test','Simulink') && exist('build_lfc_control_simulink','file')
    simResult=build_lfc_control_simulink(projectDir,p,controllers);
    for i=1:3
        [~,Ac,Ec,Cy,Cu]=closedModel(p,conditionalController(i,controllers));
        xeq=-(Ac\Ec)*p.loadStep;expected=zeros(numel(simResult(i).t),2);
        for k=1:numel(simResult(i).t)
            tau=simResult(i).t(k)-p.stepTime;
            if tau>=0,expected(k,:)=([Cy;Cu]*(xeq-expm(Ac*tau)*xeq))';end
        end
        err=max(abs(expected-[simResult(i).df(:),simResult(i).u(:)]),[],'all');
        validation.([nominal(i).name '_SimulinkMaxError'])=err;assert(err<2e-6);
    end
    validation.Simulink='PASS';
end
validation.analyticalTransferFunction='PASS';validation.ODE_LSIM='PASS';
validation.controllabilityObservability='PASS';validation.PI_Routh='PASS';validation.modifiedNyquist='PASS';
validation.saturationBound='PASS';validation.observerDemonstration='PASS';
validation.matlabVersion=version;

%% Export scientific figures and numeric frequency-response data
set(groot,'defaultAxesFontName','Arial','defaultTextFontName','Arial','defaultAxesFontSize',11);
colors=[.36 .40 .45; .04 .36 .65; .85 .32 .08];
f=figure('Visible','off','Color','w','Position',[100 100 1100 720]);tiledlayout(2,1,'TileSpacing','compact');
nexttile;hold on;for i=1:3,plot(t,nominal(i).Y(:,1),'Color',colors(i,:),'LineWidth',1.6,'DisplayName',nominal(i).name);end
 yline(.001,':','HandleVisibility','off');yline(-.001,':','HandleVisibility','off');xline(1,':','HandleVisibility','off');
 ylabel('Frequency deviation (Hz)');grid on;xlim([0 15]);legend('Location','southeast');title('A 0.01 pu load increase at t = 1 s');
nexttile;hold on;for i=2:3,plot(t,nominal(i).Y(:,2),'Color',colors(i,:),'LineWidth',1.6,'DisplayName',nominal(i).name);end
 yline(.01,'--','Steady demand','HandleVisibility','off');ylabel('Secondary command (pu)');xlabel('Time (s)');grid on;xlim([0 15]);legend('Location','southeast');saveFigure(f,out,'nominal_comparison');
f=figure('Visible','off','Color','w','Position',[100 100 1000 720]);tiledlayout(2,1,'TileSpacing','compact');
BodeData=table(w,'VariableNames',{'omega_rad_s'});
for panel=1:2,nexttile;hold on;for ci=1:2
    h=squeeze(freqresp(loops{ci},w));mag=20*log10(abs(h));ph=unwrap(angle(h))*180/pi;
    if panel==1,yy=mag;else,yy=ph;end
    semilogx(w,yy,'Color',colors(ci+1,:),'LineWidth',1.5,'DisplayName',controllers(ci).name);
    BodeData.([controllers(ci).name '_magnitude_dB'])=mag;BodeData.([controllers(ci).name '_phase_deg'])=ph;
end
set(gca,'XScale','log');xlim([.01 100]);grid on;legend('Location','southwest');xline(wc,':','2.4 rad/s','HandleVisibility','off');
if panel==1,yline(0,'--','HandleVisibility','off');ylabel('Loop magnitude (dB)');title('Secondary loop L(s) = C(s)P(s), including primary droop inside P');else,yline(-180,'--','HandleVisibility','off');ylabel('Phase (deg)');xlabel('Angular frequency (rad/s)');end
end
saveFigure(f,out,'bode_margins');writetable(BodeData,fullfile(out,'bode_data.csv'));
plot_nyquist_lfc(loops,controllers,out);
f=figure('Visible','off','Color','w','Position',[100 100 1100 500]);tiledlayout(1,2,'TileSpacing','compact');
nexttile;hold on;for i=1:3,po=nominal(i).metrics.poles;scatter(po(:,1),po(:,2),55,colors(i,:),'x','LineWidth',1.5,'DisplayName',nominal(i).name);end
xline(0,'k:','HandleVisibility','off');grid on;xlabel('Real part (1/s)');ylabel('Imaginary part (rad/s)');legend('Location','northwest');title('Closed-loop poles');
nexttile;hold on;kr=[0 logspace(-3,1.3,450)];rootsGain=rlocus(loops{1},kr);
plot(real(rootsGain.'),imag(rootsGain.'),'Color',[.55 .60 .65],'HandleVisibility','off');
po=nominal(2).metrics.poles;scatter(po(:,1),po(:,2),55,colors(2,:),'x','LineWidth',1.8,'DisplayName','Selected PI: k = 1');
% Refine the PI negative-real-axis crossing before marking marginal stability.
qa=cpi.Kp*den(1);qb=cpi.Ki*den(2)-cpi.Kp*den(3);qc=cpi.Ki*den(4);
disc=sqrt(qb*qb+4*qa*qc);if qb>=0,qcrit=2*qc/(qb+disc);else,qcrit=(-qb+disc)/(2*qa);end
gm=(den(2)*qcrit-den(4))/(p.Ksys*cpi.Kp);
critical=pole(feedback(gm*loops{1},1));
assert(min(abs(critical-1i*sqrt(qcrit)))<1e-8 && min(abs(critical+1i*sqrt(qcrit)))<1e-8);
scatter(real(critical),imag(critical),40,'r','o','DisplayName','Refined marginal gain');
xline(0,'k:','HandleVisibility','off');xlim([-16 2]);ylim([-7 7]);grid on;xlabel('Real part (1/s)');ylabel('Imaginary part (rad/s)');legend('Location','northwest');title('PI root locus: k scales both controller gains');saveFigure(f,out,'poles_root_locus');
f=figure('Visible','off','Color','w','Position',[100 100 1100 500]);tiledlayout(1,2,'TileSpacing','compact');
for ci=1:2,nexttile;hold on;fill([t;flipud(t)],[envelopeMin(:,ci);flipud(envelopeMax(:,ci))],colors(ci+1,:),'FaceAlpha',.2,'EdgeColor','none','DisplayName','43-case envelope');plot(t,nominal(ci+1).Y(:,1),'Color',colors(ci+1,:),'LineWidth',1.5,'DisplayName','Nominal');yline(.001,'k:','HandleVisibility','off');yline(-.001,'k:','HandleVisibility','off');xlim([0 15]);grid on;xlabel('Time (s)');ylabel('Frequency deviation (Hz)');title([controllers(ci).name ': fixed gains, parameter factors 0.8 to 1.2']);legend('Location','southeast');end
saveFigure(f,out,'sensitivity_envelope');
f=figure('Visible','off','Color','w','Position',[100 100 1100 700]);tiledlayout(2,1,'TileSpacing','compact');
for panel=1:2,nexttile;hold on;for n=1:4,ci=ceil(n/2);sty='--';if saturation(n).antiWindup,sty='-';end
plot(saturation(n).t,saturation(n).Y(:,panel),sty,'Color',colors(ci+1,:),'LineWidth',1.4,'DisplayName',sprintf('%s, AW=%d',saturation(n).controller,saturation(n).antiWindup));end
xline(5,':','Load on','HandleVisibility','off');xline(20,':','Load off','HandleVisibility','off');grid on;xlim([0 50]);legend('Location','northeast');if panel==1,ylabel('Frequency deviation (Hz)');title('0.03 pu load pulse; secondary command limited to +/-0.015 pu');else,ylabel('Applied command (pu)');xlabel('Time (s)');end
end;saveFigure(f,out,'saturation_antiwindup');
f=figure('Visible','off','Color','w','Position',[100 100 1100 500]);tiledlayout(1,2,'TileSpacing','compact');
for panel=1:2,nexttile;hold on;for ci=1:2,L=loops{ci};if panel==1,H=feedback(1,L);else,H=minreal(controllerTF(controllers(ci))*feedback(1,L));end
semilogx(w,20*log10(abs(squeeze(freqresp(H,w)))),'Color',colors(ci+1,:),'LineWidth',1.5,'DisplayName',controllers(ci).name);end
set(gca,'XScale','log');grid on;xlim([.01 1000]);xlabel('Angular frequency (rad/s)');ylabel('Magnitude (dB)');legend('Location','southeast');if panel==1,title('Sensitivity S = 1/(1+L)');else,title('Measurement noise to command: |C S|');end
end;saveFigure(f,out,'sensitivity_noise');
f=figure('Visible','off','Color','w','Position',[100 100 1100 500]);tiledlayout(1,2,'TileSpacing','compact');
nexttile;plot(t,[Yplace(:,1),Ylqr(:,1)],'LineWidth',1.4);grid on;xlim([0 15]);xlabel('Time (s)');ylabel('Frequency deviation (Hz)');legend('Integral pole placement','Integral LQR');title('Supplement: full-state feedback designs');
nexttile;semilogy(tobs,max(vecnorm(obsError,2,2),1e-16),'LineWidth',1.5);grid on;xlabel('Time (s)');ylabel('Estimation-error norm (mixed units)');title('Observer demonstration: no disturbance, initial mismatch');saveFigure(f,out,'state_space_extension');

summary=struct('p',p,'specifications',spec,'design',struct('wc',wc,'phaseMarginTarget',targetPM,'Cstar_real',real(Cstar),'Cstar_imag',imag(Cstar)), ...
 'controllers',controllers,'rankControllability',rankCo,'rankObservability',rankOb,'plantDenominator',den, ...
 'plantPoles',[real(eig(A)),imag(eig(A))],'PI_RouthFirstColumn',routhFirstColumn,'nominal',{[nominal.metrics]}, ...
 'robustSummary',robustSummary,'saturationSettings',satSettings,'saturationMetrics',table2struct(satTable),'courseExtension',extension,'validation',validation);
fid=fopen(fullfile(out,'research_summary.json'),'w');fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));fclose(fid);
save(fullfile(out,'research_results.mat'),'p','A','B','E','C','P','Q','controllers','loops','cls','nominal','t','robustTable','robustSummary','saturation','extension','summary');
fid=fopen(fullfile(out,'validation_summary.txt'),'w');fprintf(fid,'MATLAB %s\n',version);fprintf(fid,'%s\n',jsonencode(validation,'PrettyPrint',true));fclose(fid);
fprintf('\nRESEARCH RUN COMPLETE. Results: %s\n',out);disp(nomTable);diary off;

%% Local functions
function [A,B,E,C,P,Q]=makePlant(p)
A=[-1/p.Tg,0,-1/(p.R*p.Tg);1/p.Tt,-1/p.Tt,0;0,p.Ksys/p.Tp,-1/p.Tp];B=[1/p.Tg;0;0];E=[0;0;-p.Ksys/p.Tp];C=[0 0 1];
P=minreal(tf(ss(A,B,C,0)),1e-9);Q=minreal(tf(ss(A,E,C,0)),1e-9);
end
function C=controllerTF(c)
s=tf('s');C=c.Kp+c.Ki/s;if c.Kd>0,C=C+c.Kd*s/(1+c.Tf*s);end
end
function [sys,Ac,Ec,Cy,Cu]=closedModel(p,c)
[A,B,E,C]=makePlant(p);
if isempty(c),Ac=A;Ec=E;Cy=C;Cu=zeros(1,3);
elseif c.Kd==0,Ac=[A-B*c.Kp*C,B*c.Ki;-C,0];Ec=[E;0];Cy=[C,0];Cu=[-c.Kp*C,c.Ki];
else,D=c.Kd/c.Tf;Ac=[A-B*(c.Kp+D)*C,B*c.Ki,-B*D;-C,0,0;-C/c.Tf,0,-1/c.Tf];Ec=[E;0;0];Cy=[C,0,0];Cu=[-(c.Kp+D)*C,c.Ki,-D];end
sys=ss(Ac,Ec,[Cy;Cu],zeros(2,1));
end
function c=conditionalController(i,controllers)
if i==1,c=[];else,c=controllers(i-1);end
end
function [t,Y,X]=response(sys,p)
tpre=(0:p.dt:p.stepTime-p.dt)';tau=(0:p.dt:p.tEnd-p.stepTime)';
[ya,~,xa]=lsim(sys,p.loadStep*ones(size(tau)),tau);
t=[tpre;tau+p.stepTime];Y=[zeros(numel(tpre),size(ya,2));ya];X=[zeros(numel(tpre),size(xa,2));xa];
end
function met=measures(t,Y,p,band)
idx=t>=p.stepTime;ta=t(idx)-p.stepTime;f=Y(idx,1);u=Y(idx,2);[nadir,k]=min(f);last=find(abs(f)>band,1,'last');
rec=NaN;if isempty(last),rec=0;elseif last<numel(ta),rec=ta(last+1);end
met=struct('nadirHz',nadir,'positivePeakHz',max(f),'maxAbsHz',max(abs(f)),'nadirAfterStep_s',ta(k),'recovery_s',rec, ...
    'IAE_Hz_s',trapz(ta,abs(f)),'uPeakPu',max(abs(u)),'uEnergy_pu2_s',trapz(ta,u.^2),'finalFrequencyHz',f(end));
end
function T=metricsTable(n)
m=[n.metrics];T=table(string({n.name})',[m.nadirHz]',[m.positivePeakHz]',[m.recovery_s]',[m.IAE_Hz_s]',[m.uPeakPu]',[m.PM_deg]',[m.GM_dB]',[m.Ms]', ...
 'VariableNames',{'Controller','Nadir_Hz','PositivePeak_Hz','Recovery_s','IAE_Hz_s','PeakCommand_pu','PM_deg','GM_dB','Ms_sampled'});
end
function result=nyquistCheck(L)
[nu,de]=tfdata(L,'v');wind=zeros(1,2);
for k=1:2,ep=10^(-3-k);R=10^(3+k);path=[-1i*logspace(log10(R),log10(ep),4000),ep*exp(1i*linspace(-pi/2,pi/2,1200)),1i*logspace(log10(ep),log10(R),4000),R*exp(1i*linspace(pi/2,-pi/2,1600))];
F=1+polyval(nu,path)./polyval(de,path);phase=unwrap(angle(F));wind(k)=round((phase(end)-phase(1))/(2*pi));end
np=sum(real(pole(L))>1e-7);nz=sum(real(pole(feedback(L,1)))>1e-7);
assert(nz==np-wind(1));result=struct('openRHPPoles',np,'originPoleExcluded',true,'windingCCW',wind(1),'windingRefinedCCW',wind(2),'closedUnstablePoles',nz);
end
function [t,Y,Z]=saturationResponse(p,c,q,aw)
[A,B,E,C]=makePlant(p);t=[];Z=[];x0=zeros(5,1);edges=[0 q.on_s q.off_s q.end_s];loads=[0 q.loadPu 0];
for seg=1:3,tt=(edges(seg):q.dt:edges(seg+1))';[tt,zz]=ode45(@(~,z)rhs(z,loads(seg)),tt,x0,odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',.01));x0=zz(end,:)';if seg>1,tt(1)=[];zz(1,:)=[];end;t=[t;tt];Z=[Z;zz];end
err=-Z(:,3);ud=c.Kp*err+Z(:,4);if c.Kd>0,ud=ud+c.Kd/c.Tf*(err-Z(:,5));end;u=max(-q.limitPu,min(q.limitPu,ud));Y=[Z(:,3),u];
 function dz=rhs(z,d)
 e=-C*z(1:3);raw=c.Kp*e+z(4);dfilt=0;if c.Kd>0,raw=raw+c.Kd/c.Tf*(e-z(5));dfilt=(e-z(5))/c.Tf;end
 applied=max(-q.limitPu,min(q.limitPu,raw));iz=c.Ki*e;if aw,iz=iz+(applied-raw)/q.Taw_s;end
 dz=[A*z(1:3)+B*applied+E*d;iz;dfilt];
 end
end
function saveFigure(f,out,name)
style_lfc_figure(f);
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector');savefig(f,fullfile(out,[name '.fig']));close(f);
end
