function verification = verify_saturation_convergence()
%VERIFY_SATURATION_CONVERGENCE Refine the four saved nonlinear stress tests.
% The main study and its primary results are not modified.
% Physical equations are restated component by component for this check.

projectDir = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(projectDir, 'results');
source = jsondecode(fileread(fullfile(resultsDir, 'research_summary.json')));
p = source.p;
q = source.saturationSettings;
controllers = source.controllers;
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-13, 'MaxStep', 0.001);
% Acceptance limits are declared numerical comparison tolerances.
frequencyTolerance = 1e-6;  % Hz: one thousandth of the recovery band.
commandTolerance = 1e-6;    % pu power.
integralTolerance = 1e-6;   % pu power: the stored integral output.
recoveryTolerance = q.dt + 1e-10; % One original output sample.
rows = cell(0,1);
row = 0;

for ci = 1:numel(controllers)
    c = controllers(ci);
    for aw = [false true]
        row = row + 1;
        csvName = sprintf('saturation_%s_AW%d.csv', c.name, aw);
        original = readtable(fullfile(resultsDir, csvName));
        tTarget = original.Time_s;
        [t, z] = refineTrajectory(tTarget, p, c, q, aw, options);
        assert(numel(t) == numel(tTarget) && max(abs(t-tTarget)) < 1e-10, ...
            'The refined and original sample times are not identical.');
        [~, command] = commandValues(z, c, q);
        frequencyError = max(abs(z(:,3)-original.FrequencyDeviation_Hz));
        commandError = max(abs(command-original.Command_pu));
        integralError = max(abs(z(:,4)-original.IntegralOutput_pu));
        oldRecovery = recoveryTime(tTarget, original.FrequencyDeviation_Hz, ...
            q.off_s, source.specifications.recoveryBandHz);
        newRecovery = recoveryTime(t, z(:,3), q.off_s, ...
            source.specifications.recoveryBandHz);
        recoveryDifference = newRecovery-oldRecovery;
        passed = frequencyError <= frequencyTolerance && ...
            commandError <= commandTolerance && integralError <= integralTolerance && ...
            isfinite(oldRecovery) && isfinite(newRecovery) && ...
            abs(recoveryDifference) <= recoveryTolerance;
        rows{row} = struct('Controller', c.name, 'AntiWindup', aw, ...
            'Samples', numel(t), 'MaximumFrequencyDifference_Hz', frequencyError, ...
            'MaximumCommandDifference_pu', commandError, ...
            'MaximumIntegralOutputDifference_pu', integralError, ...
            'OriginalRecoveryAfterRemoval_s', oldRecovery, ...
            'RefinedRecoveryAfterRemoval_s', newRecovery, ...
            'RecoveryDifference_s', recoveryDifference, 'Pass', passed);
    end
end

rows = vertcat(rows{:});
verification = struct2table(rows);
reportPath = fullfile(resultsDir, 'nonlinear_convergence.txt');
fid = fopen(reportPath, 'w');
assert(fid >= 0, 'Cannot write the nonlinear convergence report.');
closer = onCleanup(@() fclose(fid));
fprintf(fid, 'NONLINEAR NUMERICAL CONVERGENCE CHECK\n');
fprintf(fid, 'Completed: %s\nMATLAB: %s\n\n', char(datetime('now')), version);
fprintf(fid, ['Scope: the four saved PI/PIDF saturation trajectories, with and without anti-windup.\n' ...
    'No parameter sweep or controller retuning was performed.\n' ...
    'The physical differential equations were restated independently, component by component.\n' ...
    'The input pulse and both segment boundaries match the main study exactly.\n\n']);
fprintf(fid, ['Original method: ode45, RelTol=1e-8, AbsTol=1e-10, MaxStep=0.01 s.\n' ...
    'Refined method: ode45, RelTol=1e-11, AbsTol=1e-13, MaxStep=0.001 s.\n' ...
    'Both methods split integration at load application and removal.\n' ...
    'Comparison uses every original CSV sample time, spaced by %.6g s.\n\n'], q.dt);
fprintf(fid, ['Declared numerical tolerances:\n' ...
    '  Maximum frequency difference <= %.3g Hz.\n' ...
    '  Maximum command difference <= %.3g pu.\n' ...
    '  Maximum integral-output difference <= %.3g pu.\n' ...
    '  Recovery-time difference <= %.6g s, one original output sample.\n\n'], ...
    frequencyTolerance, commandTolerance, integralTolerance, q.dt);
for k = 1:numel(rows)
    r = rows(k);
    status = 'FAIL'; if r.Pass, status = 'PASS'; end
    fprintf(fid, '%s, anti-windup=%d: %s\n', r.Controller, r.AntiWindup, status);
    fprintf(fid, '  Samples compared: %d\n', r.Samples);
    fprintf(fid, '  Maximum frequency difference: %.12g Hz\n', r.MaximumFrequencyDifference_Hz);
    fprintf(fid, '  Maximum command difference: %.12g pu\n', r.MaximumCommandDifference_pu);
    fprintf(fid, '  Maximum integral-output difference: %.12g pu\n', r.MaximumIntegralOutputDifference_pu);
    fprintf(fid, '  Original recovery after load removal: %.9f s\n', r.OriginalRecoveryAfterRemoval_s);
    fprintf(fid, '  Refined recovery after load removal:  %.9f s\n', r.RefinedRecoveryAfterRemoval_s);
    fprintf(fid, '  Refined minus original recovery: %.12g s\n\n', r.RecoveryDifference_s);
end
fprintf(fid, ['Recovery means remaining within +/-%.6g Hz after load removal at %.6g s,\n' ...
    'at all later sampled times through the %.6g s simulation endpoint.\n' ...
    'It is not a continuous-time crossing estimate between samples.\n\n'], ...
    source.specifications.recoveryBandHz, q.off_s, q.end_s);
if all(verification.Pass)
    fprintf(fid, 'OVERALL: PASS. All four saved trajectories agree under this numerical refinement.\n');
else
    fprintf(fid, 'OVERALL: FAIL. Inspect the failed rows before using the original trajectories.\n');
end
fprintf(fid, ['This check supports numerical convergence for these four scenarios only.\n' ...
    'It is not physical validation, a rigorous global error bound, or an independent check of every modelling assumption.\n' ...
    'The main study and its saved primary results were not modified.\n']);
disp(verification);
fprintf('Nonlinear convergence report: %s\n', reportPath);
assert(all(verification.Pass),'LFC:NonlinearConvergenceFailed', ...
    'Nonlinear convergence failed. Inspect results/nonlinear_convergence.txt before using the trajectories.');
end

function [t, z] = refineTrajectory(tTarget, p, c, q, aw, options)
edges = [0 q.on_s q.off_s q.end_s];
loads = [0 q.loadPu 0];
t = []; z = []; initial = zeros(5,1);
for seg = 1:3
    selectedTimes = tTarget(tTarget >= edges(seg) & tTarget <= edges(seg+1));
    [ts, zs] = ode45(@(~, state) physicalDynamics(state, loads(seg), p, c, q, aw), ...
        selectedTimes, initial, options);
    initial = zs(end,:)';
    if seg > 1, ts(1) = []; zs(1,:) = []; end
    t = [t; ts]; z = [z; zs]; %#ok<AGROW>
end
end

function derivative = physicalDynamics(state, loadChange, p, c, q, aw)
governor = state(1); mechanical = state(2); frequency = state(3);
error = -frequency;
[requested, applied] = commandValues(state', c, q);
integralRate = c.Ki*error;
if aw, integralRate = integralRate+(applied-requested)/q.Taw_s; end
filterRate = 0;
if c.Kd > 0, filterRate = (error-state(5))/c.Tf; end
derivative = [(applied-frequency/p.R-governor)/p.Tg; ...
    (governor-mechanical)/p.Tt; ...
    (p.Ksys*(mechanical-loadChange)-frequency)/p.Tp; ...
    integralRate; filterRate];
end

function [requested, applied] = commandValues(states, c, q)
error = -states(:,3);
requested = c.Kp*error+states(:,4);
if c.Kd > 0, requested = requested+(c.Kd/c.Tf)*(error-states(:,5)); end
applied = min(q.limitPu, max(-q.limitPu, requested));
end

function recovery = recoveryTime(t, frequency, removalTime, band)
after = t >= removalTime;
elapsed = t(after)-removalTime;
lastOutside = find(abs(frequency(after)) > band, 1, 'last');
if isempty(lastOutside)
    recovery = 0;
elseif lastOutside == numel(elapsed)
    recovery = NaN;
else
    recovery = elapsed(lastOutside+1);
end
end
