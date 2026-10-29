function summary = refine_pi_gain_margins()
%REFINE_PI_GAIN_MARGINS Refine PI gain margins without rerunning simulations.
% Uses the saved gains and 43 parameter cases. Only results/review_revision/
% is written. Existing margins, phase margins and trajectories are untouched.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'results','review_revision');if ~isfolder(out),mkdir(out);end
src=jsondecode(fileread(fullfile(root,'results','research_summary.json')));
allCases=readtable(fullfile(root,'results','parameter_sensitivity.csv'),'TextType','string');
rows=allCases(allCases.controller=="PI",:);assert(height(rows)==43);
c=src.controllers(strcmp({src.controllers.name},'PI'));p=src.p;spec=src.specifications;
assert(c.Kp>0 && c.Ki>0 && c.Kd==0);
canCheck=~isempty(which('allmargin')) && ~isempty(which('tf'));
records=cell(height(rows),1);checks=cell(0,1);factorNames={'Tg','Tt','Tp','Ksys','R'};
for i=1:height(rows)
    pv=p;
    for j=1:numel(factorNames)
        name=factorNames{j};pv.(name)=p.(name)*rows.([name '_factor'])(i);
    end
    a3=pv.Tg*pv.Tt*pv.Tp;
    a2=pv.Tg*pv.Tt+pv.Tg*pv.Tp+pv.Tt*pv.Tp;
    a1=pv.Tg+pv.Tt+pv.Tp;a0=1+pv.Ksys/pv.R;
    hurwitzSlack=a2*a1-a3*a0;assert(hurwitzSlack>0);
    % Im L(jw)=0 reduces to A*q^2+b*q-c0=0, q=w^2.
    % A and c0 are positive, giving exactly one positive real root.
    qa=c.Kp*a3;qb=c.Ki*a2-c.Kp*a1;qc=c.Ki*a0;
    disc=sqrt(qb*qb+4*qa*qc);
    if qb>=0,q=2*qc/(qb+disc);else,q=(-qb+disc)/(2*qa);end
    w=sqrt(q);gm=(a2*q-a0)/(pv.Ksys*c.Kp);gmDb=20*log10(gm);
    den=[a3 a2 a1 a0];z=1i*w;
    Ljw=pv.Ksys*(c.Kp+c.Ki/z)/polyval(den,z);
    phaseResidual=abs(imag(Ljw));criticalResidual=abs(1+gm*Ljw);
    assert(gm>0 && real(Ljw)<0 && phaseResidual<1e-12 && criticalResidual<1e-12);
    % A gain-scaled PI quartic must contain the predicted imaginary pair.
    charpoly=[a3 a2 a1 a0+gm*pv.Ksys*c.Kp gm*pv.Ksys*c.Ki];
    criticalPoles=roots(charpoly);
    poleResidual=max(min(abs(criticalPoles-1i*w)),min(abs(criticalPoles+1i*w)));
    assert(poleResidual<1e-9);
    legacyGm=rows.GM_dB(i);
    common=logical(rows.stable(i)) && rows.recovery_s(i)<=spec.recoveryTime_s ...
        && rows.maxAbsHz(i)<=spec.maxAbsoluteDeviationHz ...
        && rows.uPeakPu(i)<=spec.maxSecondaryCommandPu ...
        && rows.PM_deg(i)>=spec.minPhaseMarginDeg;
    legacyPass=legacyGm>=spec.minGainMarginDb;refinedPass=gmDb>=spec.minGainMarginDb;
    records{i}=struct('CaseID',rows.caseID(i),'Tg_factor',rows.Tg_factor(i), ...
        'Tt_factor',rows.Tt_factor(i),'Tp_factor',rows.Tp_factor(i), ...
        'Ksys_factor',rows.Ksys_factor(i),'R_factor',rows.R_factor(i), ...
        'LegacyGM_dB',legacyGm,'RefinedGM_dB',gmDb,'RefinedMinusLegacy_dB',gmDb-legacyGm, ...
        'RefinedGM_linear',gm,'RefinedPhaseCrossing_rad_s',w, ...
        'LegacyGMTargetPass',legacyPass,'RefinedGMTargetPass',refinedPass, ...
        'LegacyAllTargetsPass',common&&legacyPass,'RefinedAllTargetsPass',common&&refinedPass, ...
        'PrimaryPlantHurwitzSlack',hurwitzSlack, ...
        'PhaseCrossingImaginaryResidual',phaseResidual,'CriticalLoopResidual',criticalResidual, ...
        'CriticalPolePairResidual_per_s',poleResidual);
    if canCheck
        L=tf(pv.Ksys*[c.Kp c.Ki],conv([1 0],den));
        [mg,~,mw]=margin(L);am=allmargin(L);
        valid=find(isfinite(am.GainMargin) & am.GainMargin>0);
        assert(numel(valid)==1,'Unexpected number of finite positive gain margins.');
        ag=am.GainMargin(valid);aw=am.GMFrequency(valid);
        atMatlab=pv.Ksys*(c.Kp+c.Ki/(1i*aw))/polyval(den,1i*aw);
        checks{end+1}=struct('CaseID',rows.caseID(i),'MarginGM_dB',20*log10(mg), ...
            'MarginPhaseCrossing_rad_s',mw,'AllmarginGM_dB',20*log10(ag), ...
            'AllmarginPhaseCrossing_rad_s',aw,'AnalyticGM_dB',gmDb, ...
            'AllmarginMinusAnalytic_dB',20*log10(ag)-gmDb, ...
            'ImaginaryLoopAtAllmarginCrossing',imag(atMatlab)); %#ok<AGROW>
    end
end
tableOut=struct2table(vertcat(records{:}));
writetable(tableOut,fullfile(out,'pi_gain_margin_refinement.csv'));
[~,worstIndex]=min(tableOut.RefinedGM_dB);
[~,largestIndex]=max(abs(tableOut.RefinedMinusLegacy_dB));
changed=tableOut.CaseID(tableOut.LegacyGMTargetPass~=tableOut.RefinedGMTargetPass ...
    | tableOut.LegacyAllTargetsPass~=tableOut.RefinedAllTargetsPass);
summary=struct('executionEngine','MATLAB','matlabVersion',version,'inputSource','Saved research_summary.json and parameter_sensitivity.csv', ...
    'controller',c,'caseCount',height(rows),'gainMarginTarget_dB',spec.minGainMarginDb, ...
    'nominal',records{find(rows.caseID==1,1)},'worst',records{worstIndex}, ...
    'largestAbsoluteCorrection',records{largestIndex}, ...
    'legacyGMTargetPassCount',sum(tableOut.LegacyGMTargetPass), ...
    'refinedGMTargetPassCount',sum(tableOut.RefinedGMTargetPass), ...
    'legacyAllTargetsPassCount',sum(tableOut.LegacyAllTargetsPass), ...
    'refinedAllTargetsPassCount',sum(tableOut.RefinedAllTargetsPass), ...
    'changedClassificationCaseIDs',changed','phaseMarginsAndTrajectoriesChanged',false, ...
    'minimumPrimaryPlantHurwitzSlack',min(tableOut.PrimaryPlantHurwitzSlack), ...
    'maxNegativeAxisImaginaryResidual',max(tableOut.PhaseCrossingImaginaryResidual), ...
    'maxCriticalLoopResidual',max(tableOut.CriticalLoopResidual), ...
    'analyticMethod','q=w^2; Kp*a3*q^2+(Ki*a2-Kp*a1)*q-Ki*a0=0; GM=(a2*q-a0)/(Ksys*Kp)', ...
    'analyticScope','PI controller with positive gains and the retained third-order primary-droop plant; exact algebra evaluated in floating point', ...
    'legacyProvenance','run_research.m parameter sweep calls margin(Pv*controllerTF(...)); its separate plotted frequency grid is not passed to margin', ...
    'allmarginCrossCheckExecuted',canCheck,'cases',vertcat(records{:}));
if canCheck
    summary.allmarginCrossChecks=vertcat(checks{:});
    writetable(struct2table(summary.allmarginCrossChecks),fullfile(out,'pi_gain_margin_matlab_comparison.csv'));
else
    summary.allmarginCrossCheckStatus='Control System Toolbox unavailable; closed-form refinement completed without a fresh margin/allmargin comparison';
end
fid=fopen(fullfile(out,'pi_gain_margin_summary.json'),'w');assert(fid>=0);closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));clear closeFile
fid=fopen(fullfile(out,'pi_gain_margin_findings.txt'),'w');assert(fid>=0);closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'PI GAIN-MARGIN PRECISION CHECK\n\n');
fprintf(fid,'The saved gains, physical model, phase margins and time responses are retained.\n');
fprintf(fid,'The 43 PI cases are refined analytically without rerunning any trajectory.\n');
fprintf(fid,'Nominal: legacy %.12f dB; refined %.12f dB; phase crossing %.12f rad/s.\n', ...
    summary.nominal.LegacyGM_dB,summary.nominal.RefinedGM_dB,summary.nominal.RefinedPhaseCrossing_rad_s);
fprintf(fid,'Worst case %d: legacy %.12f dB; refined %.12f dB; phase crossing %.12f rad/s.\n', ...
    summary.worst.CaseID,summary.worst.LegacyGM_dB,summary.worst.RefinedGM_dB,summary.worst.RefinedPhaseCrossing_rad_s);
fprintf(fid,'Largest correction, case %d: %.12f dB.\n',summary.largestAbsoluteCorrection.CaseID,summary.largestAbsoluteCorrection.RefinedMinusLegacy_dB);
fprintf(fid,'GM target counts: %d/43 before and %d/43 after.\n',summary.legacyGMTargetPassCount,summary.refinedGMTargetPassCount);
fprintf(fid,'All-target counts: %d/43 before and %d/43 after.\n',summary.legacyAllTargetsPassCount,summary.refinedAllTargetsPassCount);
if isempty(changed),changedText='none';else,changedText=mat2str(changed');end
fprintf(fid,'Changed classification case IDs: %s.\n',changedText);
fprintf(fid,'The all-target check retains saved stability, recovery, maximum deviation, command and phase-margin classifications.\n');
fprintf(fid,'This refinement changes gain-margin precision, not the model, controller tuning or physical evidence.\n');
fprintf(fid,'Fresh MATLAB margin/allmargin comparison executed: %d.\n',canCheck);
if canCheck
    k=find([summary.allmarginCrossChecks.CaseID]==summary.worst.CaseID,1);
    worstCheck=summary.allmarginCrossChecks(k);
    fprintf(fid,'Fresh worst-case allmargin: %.12f dB at %.12f rad/s.\n', ...
        worstCheck.AllmarginGM_dB,worstCheck.AllmarginPhaseCrossing_rad_s);
    fprintf(fid,'Its returned crossover has directly evaluated imaginary-loop residual %.12g.\n', ...
        worstCheck.ImaginaryLoopAtAllmarginCrossing);
    fprintf(fid,'The analytic worst-case crossing has imaginary-loop residual %.12g.\n',summary.worst.PhaseCrossingImaginaryResidual);
end
fprintf(fid,'The comparison CSV reports returned values and directly evaluated crossing residuals; it does not infer an internal software algorithm.\n');
fprintf(fid,'The legacy sweep uses margin() directly. Its separate plotted frequency grid is not the source supplied to that call.\n');
clear closeFile
assert(isempty(changed),'A PI classification changed; inspect the new comparison before updating narrative results.');
fprintf('PI margin refinement complete: nominal %.12f dB, worst %.12f dB, classifications unchanged.\n', ...
    summary.nominal.RefinedGM_dB,summary.worst.RefinedGM_dB);
end
