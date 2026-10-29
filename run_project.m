%% Reproduce the complete ENGN6223 load frequency control study
% Run from any working folder. Generated results in this project are replaced.
% Additional diagnostics are implemented in matlab/.
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot,'matlab'));
run(fullfile(projectRoot,'matlab','run_research.m'));
compare_parameter_recovery();
verify_saturation_convergence();
analyse_recovery_exceptions();
analyse_phase_margin();
refine_pi_gain_margins();
analyse_design_sensitivity();
clarify_design_horizon();
if strcmp(validation.Simulink,'PASS')
    export_simulink_diagrams(projectRoot);
end
fprintf('\nProject reproduction finished. Read results/validation_summary.txt, results/deepening/ and results/review_revision/.\n');
