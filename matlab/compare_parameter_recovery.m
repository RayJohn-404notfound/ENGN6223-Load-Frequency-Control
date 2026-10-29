function comparison=compare_parameter_recovery()
%COMPARE_PARAMETER_RECOVERY Compare controllers at matching parameter cases.
root=fileparts(fileparts(mfilename('fullpath')));
T=readtable(fullfile(root,'results','parameter_sensitivity.csv'));
piRows=T(strcmp(T.controller,'PI'),:);
pidRows=T(strcmp(T.controller,'PIDF'),:);
assert(isequal(piRows.caseID,pidRows.caseID),'Parameter case IDs do not match.');
fasterPI=piRows.recovery_s<pidRows.recovery_s;
fasterPID=pidRows.recovery_s<piRows.recovery_s;
comparison=struct('basis','Pairwise comparison at matching parameter case IDs.', ...
    'cases',height(piRows),'PIDF_faster',sum(fasterPID),'PI_faster',sum(fasterPI), ...
    'equalRecovery',sum(piRows.recovery_s==pidRows.recovery_s));
exceptions=table(piRows.caseID(fasterPI),piRows.recovery_s(fasterPI), ...
    pidRows.recovery_s(fasterPI),'VariableNames',{'caseID','PI_s','PIDF_s'});
comparison.PI_faster_cases=table2struct(exceptions);
fid=fopen(fullfile(root,'results','recovery_comparison.json'),'w');
assert(fid~=-1,'Could not open recovery comparison output.');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(comparison,'PrettyPrint',true));
fprintf('PIDF recovers faster in %d/%d matched parameter cases.\n', ...
    comparison.PIDF_faster,comparison.cases);
end
