function summary=summarise_design_sensitivity()
%SUMMARISE_DESIGN_SENSITIVITY Rebuild compact summaries without new responses.
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'results','design_sensitivity');
summary=jsondecode(fileread(fullfile(out,'design_sensitivity_summary.json')));
A=readtable(fullfile(out,'design_summary.csv'),'TextType','string');T=readtable(fullfile(out,'all_plant_cases.csv'),'TextType','string');
P=readtable(fullfile(out,'matched_comparisons.csv'),'TextType','string');N=readtable(fullfile(out,'nominal_metrics.csv'),'TextType','string');
% A censored case cannot be omitted from a quantity called worst recovery.
A.WorstRecovery_s(A.HorizonLimitedCases>0)=NaN;writetable(A,fullfile(out,'design_summary.csv'));
summary.aggregate=table2struct(A);summary.multipleUnityCrossingCases=sum(T.UnityCrossingCount>1);
summary.singleUnityCrossingCases=sum(T.UnityCrossingCount==1);
summary.horizonBoundary='WorstRecovery_s is null/NaN if any original-horizon response is censored. Separate continuation never replaces the primary 30 s metric.';
E=P(logical(P.Eligible),:);summary.pairedNominalPIDFFaster=sum(E.PIDF_NominalRecovery_s<E.PI_NominalRecovery_s);
summary.pairedNominalPIDFLowerIAE=sum(E.PIDF_NominalIAE_Hz_s<E.PI_NominalIAE_Hz_s);
summary.pairedNominalPIDFHigherNoise=sum(E.PIDF_Noise100>E.PI_Noise100);
summary.pairingBoundary='15 eligible design pairs reuse five PI controllers across three N values; their 645 plant comparisons are repeated deterministic comparisons, not independent samples.';
keys=unique(E(:,{'Wc_rad_s','TargetPM_deg'}),'rows');rows=cell(height(keys),1);
for i=1:height(keys)
 r=E(E.Wc_rad_s==keys.Wc_rad_s(i) & E.TargetPM_deg==keys.TargetPM_deg(i),:);assert(height(r)==3);
 rows{i}=struct('Wc_rad_s',keys.Wc_rad_s(i),'TargetPM_deg',keys.TargetPM_deg(i),'PI_Recovery_s',r.PI_NominalRecovery_s(1), ...
  'PIDF_RecoveryMin_s',min(r.PIDF_NominalRecovery_s),'PIDF_RecoveryMax_s',max(r.PIDF_NominalRecovery_s), ...
  'PI_Noise100',r.PI_Noise100(1),'PIDF_Noise100Min',min(r.PIDF_Noise100),'PIDF_Noise100Max',max(r.PIDF_Noise100), ...
  'PI_AllTargetCount',r.PI_AllPassCases(1),'PIDF_AllTargetCountMin',min(r.PIDF_AllPassCases),'PIDF_AllTargetCountMax',max(r.PIDF_AllPassCases), ...
  'PI_IAE_Hz_s',r.PI_NominalIAE_Hz_s(1),'PIDF_IAEMin_Hz_s',min(r.PIDF_NominalIAE_Hz_s),'PIDF_IAEMax_Hz_s',max(r.PIDF_NominalIAE_Hz_s));
end
K=struct2table(vertcat(rows{:}));writetable(K,fullfile(out,'matched_target_summary.csv'));summary.compactMatchedTargets=table2struct(K);
F=N(N.Type=="PIDF",:);keys=unique(F(:,{'Wc_rad_s','TargetPM_deg'}),'rows');rows=cell(height(keys),1);
for i=1:height(keys)
 r=sortrows(F(F.Wc_rad_s==keys.Wc_rad_s(i) & F.TargetPM_deg==keys.TargetPM_deg(i),:),'FilterN');assert(isequal(r.FilterN,[5;10;20]));
 rows{i}=struct('Wc_rad_s',keys.Wc_rad_s(i),'TargetPM_deg',keys.TargetPM_deg(i), ...
  'RecoveryN5_s',r.Recovery_s(1),'RecoveryN10_s',r.Recovery_s(2),'RecoveryN20_s',r.Recovery_s(3), ...
  'IAEN5_Hz_s',r.IAE_Hz_s(1),'IAEN10_Hz_s',r.IAE_Hz_s(2),'IAEN20_Hz_s',r.IAE_Hz_s(3), ...
  'NoiseN5',r.CommandNoiseGain100(1),'NoiseN10',r.CommandNoiseGain100(2),'NoiseN20',r.CommandNoiseGain100(3), ...
  'MsN5',r.Ms(1),'MsN10',r.Ms(2),'MsN20',r.Ms(3));
end
F=struct2table(vertcat(rows{:}));writetable(F,fullfile(out,'filter_redesign_summary.csv'));summary.filterRedesigns=table2struct(F);
f=fopen(fullfile(out,'design_sensitivity_summary.json'),'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(summary,'PrettyPrint',true));fclose(f);
f=fopen(fullfile(out,'interpretation.txt'),'w');assert(f>=0);clean=onCleanup(@()fclose(f));
fprintf(f,'BOUNDED INTERPRETATION OF THE DESIGN GRID\n\n');
fprintf(f,'All 36 planned candidates remain registered. Eight PI and 18 PIDF candidates have positive gains and stable nominal closed loops.\n');
fprintf(f,'The exclusions arise from the specified tuning closure. They do not show that general PIDF cannot realise the excluded target pairs.\n');
fprintf(f,'Fifteen eligible matched pairs occupy five shared crossover/margin cells. The same PI design is reused for each of three PIDF filter settings.\n');
fprintf(f,'PIDF gives shorter nominal recovery, lower IAE and greater command-noise gain in all 15 pairs. The 55-degree cell improves recovery by only 7–8 ms.\n');
fprintf(f,'Forty-five controller-plant combinations have three unity-gain crossings. Margins use every crossing, and the nearest signed margin governs screening.\n');
fprintf(f,'All 1,118 evaluated plant/controller combinations are stable. One original-horizon recovery is censored; its worst-recovery summary remains NaN.\n');
fprintf(f,'At the retained 2.4 rad/s, 65-degree target, N5/10/20 produces recovery 2.972/2.991/3.002 s and command-noise gain 0.673/1.100/1.846.\n');
fprintf(f,'The same redesigns give IAE 0.020094/0.019879/0.019799 Hz s. Lower noise and faster band recovery therefore do not minimise every metric simultaneously.\n');
fprintf(f,'Changing N re-solves all PIDF gains. This comparison cannot isolate the filter alone.\n');
fprintf(f,'The original design remains the benchmark. This finite exploratory grid does not establish an optimum or physical performance.\n');
end
