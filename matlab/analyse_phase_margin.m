function summary=analyse_phase_margin()
%ANALYSE_PHASE_MARGIN Explain matched phase margins and separate gain effects.
% Adds evidence without retuning or overwriting the original controller study.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'results','deepening');if ~isfolder(out),mkdir(out);end
src=jsondecode(fileread(fullfile(root,'results','research_summary.json')));
p=src.p;cs=src.controllers;wc=src.design.wc;
[A,B,E,C,P,Q]=plant(p); %#ok<ASGLU>
tau=(0:.001:60)';names={'PI','PIDF'};
rows=cell(2,1);modalRows={};k=0;responses=zeros(numel(tau),2);
sens=zeros(3000,2);omega=logspace(-2,2,3000)';loops=cell(2,1);
for i=1:2
    c=cs(i);[Ac,Ec,Cy,Cu]=closed(A,B,E,C,c);
    [V,D]=eig(Ac);lambda=diag(D);xeq=-(Ac\Ec)*p.loadStep;
    coeff=-(Cy*V).'.*(V\xeq);yInf=Cy*xeq;
    ym=real(exp(tau*lambda.')*coeff+yInf);
    [y,~,~]=lsim(ss(Ac,Ec,[Cy;Cu],zeros(2,1)),p.loadStep*ones(size(tau)),tau);
    assert(max(abs(ym-y(:,1)))<1e-10,'Modal and state responses disagree.');
    responses(:,i)=y(:,1);
    Ci=controller(c);L=minreal(P*Ci,1e-9);loops{i}=L;S=feedback(1,L);
    [ms,wp]=getPeakGain(S,1e-8);
    % Independent local refinement around the reported peak.
    wr=exp(fminbnd(@(v)-abs(evalfr(S,1i*exp(v))),log(wp/1.5),log(wp*1.5), ...
        optimset('TolX',1e-12,'Display','off')));
    mr=abs(evalfr(S,1i*wr));assert(abs(ms-mr)<1e-6);
    am=allmargin(L);Lc=evalfr(L,1i*wc);Lp=evalfr(L,1i*wp);
    dist=abs(1+Lp);sc=abs(evalfr(S,1i*wc));
    assert(abs(ms*dist-1)<1e-6);
    assert(abs(sc-1/(2*sind(src.design.phaseMarginTarget/2)))<1e-10);
    sens(:,i)=abs(squeeze(freqresp(S,omega)));
    pair=find(imag(lambda)>1e-6);reals=find(abs(imag(lambda))<1e-6);
    pairSignal=real(exp(tau*lambda(pair).')*coeff(pair))*2;
    realSignal=real(exp(tau*lambda(reals).')*coeff(reals));
    [rec,recRoot]=recovery(tau,ym,.001,lambda,coeff,yInf);
    rows{i}=struct('controller',names{i},'PM_deg',am.PhaseMargin, ...
        'crossover_rad_s',am.PMFrequency,'Ms_peak',ms,'Ms_peak_rad_s',wp, ...
        'Ms_fminbnd',mr,'Ms_refined_rad_s',wr, ...
        'sensitivityAtCrossover',sc,'distanceAtCrossover',abs(1+Lc), ...
        'minimumNyquistDistance',dist,'closestLoopReal',real(Lp),'closestLoopImag',imag(Lp), ...
        'commonLoopReal',real(Lc),'commonLoopImag',imag(Lc), ...
        'recoverySampled_s',rec,'recoveryRoot_s',recRoot, ...
        'dominantDecayRate_s_inv',max(real(lambda)), ...
        'dominantPairReal',real(lambda(pair)),'dominantPairImag',imag(lambda(pair)), ...
        'pairCoefficientAbs_Hz',abs(coeff(pair)), ...
        'pairCoefficientPhase_rad',angle(coeff(pair)), ...
        'modalStateMaxDifference_Hz',max(abs(ym-y(:,1))));
    for j=1:numel(lambda)
        k=k+1;modalRows{k}=struct('controller',names{i},'poleReal_s_inv',real(lambda(j)), ...
            'poleImag_rad_s',imag(lambda(j)),'coefficientReal_Hz',real(coeff(j)), ...
            'coefficientImag_Hz',imag(coeff(j)),'coefficientMagnitude_Hz',abs(coeff(j)), ...
            'coefficientPhase_rad',angle(coeff(j)));
    end
    writetable(table(tau,ym,pairSignal,realSignal,'VariableNames', ...
        {'TimeAfterStep_s','Frequency_Hz','OscillatoryPair_Hz','RealModes_Hz'}), ...
        fullfile(out,[names{i} '_nominal_modes.csv']));
end
nominal=vertcat(rows{:});
writetable(struct2table(vertcat(modalRows{:})),fullfile(out,'nominal_modal_coefficients.csv'));
writetable(table(omega,sens(:,1),sens(:,2),'VariableNames', ...
    {'omega_rad_s','PI_S_magnitude','PIDF_S_magnitude'}),fullfile(out,'sensitivity_peak_data.csv'));

%% Full 2x2x2 controlled gain comparison; no equal-margin claim for these variants
% The three factors select the published PI or PIDF proportional/integral gain,
% and zero or the published filtered derivative term. Tf stays fixed when used.
bits=dec2bin(0:7,3)-'0';variants=cell(8,1);vResponses=zeros(numel(tau),8);
for i=1:8
    c=cs(1);c.Kp=cs(1+bits(i,1)).Kp;c.Ki=cs(1+bits(i,2)).Ki;
    c.Kd=bits(i,3)*cs(2).Kd;c.Tf=cs(2).Tf;
    [Ac,Ec,Cy,Cu]=closed(A,B,E,C,c);[V,D]=eig(Ac);lambda=diag(D);
    xeq=-(Ac\Ec)*p.loadStep;coeff=-(Cy*V).'.*(V\xeq);yInf=Cy*xeq;
    y=real(exp(tau*lambda.')*coeff+yInf);
    xx=real(xeq.'-(exp(tau*lambda.').*(V\xeq).')*V.');u=xx*Cu.';
    [rec,recRoot]=recovery(tau,y,.001,lambda,coeff,yInf);
    fineTau=(0:.0001:60)';
    fineY=real(exp(fineTau*lambda.')*coeff+yInf);
    [~,fineRoot]=recovery(fineTau,fineY,.001,lambda,coeff,yInf);
    assert(abs(fineRoot-recRoot)<1e-8,'Recovery is sensitive to the sampling grid.');
    Ci=controller(c);L=P*Ci;am=allmargin(L);[ms,wp]=getPeakGain(feedback(1,L),1e-8);
    assert(all(real(lambda)<0),'A diagnostic gain variant is unstable.');
    variants{i}=struct('variantID',i,'P_uses_PIDF_gain',logical(bits(i,1)), ...
        'I_uses_PIDF_gain',logical(bits(i,2)),'D_enabled',logical(bits(i,3)), ...
        'Kp',c.Kp,'Ki',c.Ki,'Kd',c.Kd,'Tf_s',c.Tf, ...
        'stable',true,'recovery_s',rec,'recoveryRoot_s',recRoot, ...
        'recoveryRootFine_s',fineRoot,'recoveryRootGridDifference_s',abs(fineRoot-recRoot), ...
        'nadir_Hz',min(y),'positivePeak_Hz',max(y),'IAE_29s_Hz_s',trapz(tau(tau<=29),abs(y(tau<=29))), ...
        'uPeak_pu',max(abs(u)),'PM_deg',min(am.PhaseMargin), ...
        'crossover_rad_s',am.PMFrequency,'Ms',ms,'Ms_peak_rad_s',wp, ...
        'noiseCommandGain100_pu_Hz',abs(evalfr(Ci*feedback(1,L),100i)), ...
        'dominantDecayRate_s_inv',max(real(lambda)));
    vResponses(:,i)=y;
end
v=vertcat(variants{:});writetable(struct2table(v),fullfile(out,'gain_ablation.csv'));
assert(abs(v(1).recovery_s-src.nominal(2).recovery_s)<1e-9);
assert(abs(v(8).recovery_s-src.nominal(3).recovery_s)<1e-9);
% Factorial effect is a descriptive finite difference, not a universal causal coefficient.
metrics={'recoveryRoot_s','IAE_29s_Hz_s','Ms','noiseCommandGain100_pu_Hz'};
effects=struct;
for j=1:numel(metrics)
    vals=[v.(metrics{j})]';
    effects.(metrics{j})=struct('D_at_PI_PI',vals(2)-vals(1), ...
        'I_at_PI_noD',vals(3)-vals(1),'D_at_PI_highI',vals(4)-vals(3), ...
        'I_at_PI_withD',vals(4)-vals(2),'I_D_interaction_at_PI_Kp',vals(4)-vals(3)-vals(2)+vals(1), ...
        'P_at_highI_withD',vals(8)-vals(4));
end
summary=struct('nominal',nominal,'ablation',v,'finiteDifferences',effects, ...
    'sharedCrossoverFormula','|S(jwc)| = 1 / [2 sin(PM/2)] when |L(jwc)| = 1', ...
    'ablationBoundary','Diagnostic fixed-gain substitutions; margins and crossover are recalculated, not constrained to match.', ...
    'timeConvention','All response times are seconds after the load step; load magnitude 0.01 pu.', ...
    'peakMethod','getPeakGain(S,1e-8), checked with independent local fminbnd refinement.', ...
    'recoveryValidation','Last crossings bracketed on 0.001 s and 0.0001 s grids, then refined by fzero. Numerical convergence is not a formal proof of crossing completeness.', ...
    'noiseInjection','Additive measurement noise in the secondary frequency feedback only; the reported gain is at 100 rad/s, not an integrated stochastic-noise metric.', ...
    'matlabVersion',version);
fid=fopen(fullfile(out,'phase_margin_summary.json'),'w');fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));fclose(fid);
writetable(array2table([tau,vResponses],'VariableNames', ...
    [{'TimeAfterStep_s'},arrayfun(@(x)sprintf('Variant%d_Hz',x),1:8,'UniformOutput',false)]), ...
    fullfile(out,'gain_ablation_responses.csv'));

%% Scientific figures
set(groot,'defaultAxesFontName','Arial','defaultTextFontName','Arial','defaultAxesFontSize',11);
colors=[.06 .36 .67;.87 .30 .08];
f=figure('Visible','off','Color','w','Position',[100 100 1200 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;
for i=1:2
    semilogx(omega,sens(:,i),'Color',colors(i,:),'LineWidth',1.7,'DisplayName',names{i});
    plot(nominal(i).Ms_peak_rad_s,nominal(i).Ms_peak,'o','Color',colors(i,:), ...
        'MarkerFaceColor',colors(i,:),'HandleVisibility','off');
end
set(gca,'XScale','log');xline(wc,'k:','Shared crossover','HandleVisibility','off');
plot(wc,nominal(1).sensitivityAtCrossover,'ks','MarkerFaceColor','k','DisplayName','Shared value at crossover');
xlim([.3 30]);ylim([0 2.1]);grid on;xlabel('Angular frequency (rad/s)');ylabel('|S(j\omega)|');
title('Equal value at crossover; different peaks');legend('Location','southwest');
nexttile;hold on;
for i=1:2
    z=squeeze(freqresp(loops{i},logspace(-1,2,3000)));
    plot(real(z),imag(z),'Color',colors(i,:),'LineWidth',1.6,'DisplayName',names{i});
    zmin=nominal(i).closestLoopReal+1i*nominal(i).closestLoopImag;
    plot([-1 real(zmin)],[0 imag(zmin)],'--','Color',colors(i,:),'HandleVisibility','off');
    plot(real(zmin),imag(zmin),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:),'HandleVisibility','off');
end
plot(-1,0,'rx','MarkerSize',10,'LineWidth',2,'DisplayName','Critical point -1');
plot(nominal(1).commonLoopReal,nominal(1).commonLoopImag,'ks','MarkerFaceColor','k','DisplayName','Shared crossover');
axis equal;xlim([-1.2 .2]);ylim([-1.2 .3]);grid on;
xlabel('Real L(j\omega)');ylabel('Imaginary L(j\omega)');
title('Closest approach differs away from crossover');legend('Location','northwest');
saveFigure(f,out,'phase_margin_mechanism');

f=figure('Visible','off','Color','w','Position',[100 100 1200 470]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for i=1:2
    T=readtable(fullfile(out,[names{i} '_nominal_modes.csv']));
    nexttile;hold on;plot(T.TimeAfterStep_s,T.Frequency_Hz,'k','LineWidth',1.8,'DisplayName','Total frequency');
    plot(T.TimeAfterStep_s,T.OscillatoryPair_Hz,'Color',colors(i,:),'LineWidth',1.4,'DisplayName','Complex pole pair');
    plot(T.TimeAfterStep_s,T.RealModes_Hz,'--','Color',[.45 .45 .45],'LineWidth',1.3,'DisplayName','Real modes');
    yline(.001,'k:','HandleVisibility','off');yline(-.001,'k:','HandleVisibility','off');
    xlim([0 8]);ylim([-.04 .04]);grid on;xlabel('Time after load step (s)');ylabel('Frequency contribution (Hz)');
    title([names{i} ': poles and coefficients both matter']);legend('Location','northeast');
end
saveFigure(f,out,'nominal_modal_decomposition');

f=figure('Visible','off','Color','w','Position',[100 100 1200 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
labels={'Original PI','D only','I only','I and D'};co=lines(4);
nexttile;hold on;
for i=1:4,plot(tau,vResponses(:,i),'LineWidth',1.4,'Color',co(i,:),'DisplayName',labels{i});end
xlim([0 10]);grid on;yline(.001,'k:','HandleVisibility','off');yline(-.001,'k:','HandleVisibility','off');
xlabel('Time after load step (s)');ylabel('Frequency deviation (Hz)');
title('Fixed Kp: separate integral and derivative changes');legend('Location','southeast');
nexttile;bar(categorical(labels,labels),[v(1:4).recoveryRoot_s]);ylabel('Recovery into +/-0.001 Hz (s)');grid on;
title('Substitution experiment; margins are not matched');
saveFigure(f,out,'gain_ablation');
fprintf('Phase-margin deepening complete.\n');disp(struct2table(nominal));disp(struct2table(v));
end

function [A,B,E,C,P,Q]=plant(p)
A=[-1/p.Tg,0,-1/(p.R*p.Tg);1/p.Tt,-1/p.Tt,0;0,p.Ksys/p.Tp,-1/p.Tp];
B=[1/p.Tg;0;0];E=[0;0;-p.Ksys/p.Tp];C=[0 0 1];
P=minreal(tf(ss(A,B,C,0)),1e-9);Q=minreal(tf(ss(A,E,C,0)),1e-9);
end
function ci=controller(c)
s=tf('s');ci=c.Kp+c.Ki/s;if c.Kd>0,ci=ci+c.Kd*s/(1+c.Tf*s);end
end
function [Ac,Ec,Cy,Cu]=closed(A,B,E,C,c)
if c.Kd==0,Ac=[A-B*c.Kp*C B*c.Ki;-C 0];Ec=[E;0];Cy=[C 0];Cu=[-c.Kp*C c.Ki];
else,D=c.Kd/c.Tf;Ac=[A-B*(c.Kp+D)*C B*c.Ki -B*D;-C 0 0;-C/c.Tf 0 -1/c.Tf];Ec=[E;0;0];Cy=[C 0 0];Cu=[-(c.Kp+D)*C c.Ki -D];end
end
function [sampled,refined]=recovery(t,y,band,lambda,coeff,constant)
last=find(abs(y)>band,1,'last');
if isempty(last),sampled=0;refined=0;return;end
assert(last<numel(t),'Horizon does not contain recovery.');sampled=t(last+1);
fun=@(v)abs(real(sum(coeff.*exp(lambda*v))+constant))-band;
refined=fzero(fun,[t(last),t(last+1)]);
% This decreasing absolute modal envelope excludes re-entry after t(end).
assert(sum(abs(coeff).*exp(real(lambda)*t(end)))+abs(constant)<band/100);
end
function saveFigure(f,out,name)
style_lfc_figure(f);
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector');savefig(f,fullfile(out,[name '.fig']));close(f);
end
