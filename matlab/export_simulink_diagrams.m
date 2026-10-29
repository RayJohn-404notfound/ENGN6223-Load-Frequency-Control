function export_simulink_diagrams(projectRoot)
%EXPORT_SIMULINK_DIAGRAMS Export the saved nominal comparison diagram.
if nargin==0,projectRoot=fileparts(fileparts(mfilename('fullpath')));end
model='lfc_control_comparison';
load_system(fullfile(projectRoot,'model',[model '.slx']));
out=fullfile(projectRoot,'results');
print(['-s' model],'-dpng','-r150',fullfile(out,'simulink_overview.png'));
branches=find_system(model,'SearchDepth',1,'BlockType','SubSystem');
for k=1:numel(branches)
    print(['-s' branches{k}],'-dpng','-r150', ...
        fullfile(out,sprintf('simulink_branch_%d.png',k)));
end
end
