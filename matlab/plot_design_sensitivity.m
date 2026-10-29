function plot_design_sensitivity()
%PLOT_DESIGN_SENSITIVITY Rebuild scientific figure from saved grid results.
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'results','design_sensitivity');
D=readtable(fullfile(out,'design_register.csv'));T=readtable(fullfile(out,'nominal_metrics.csv'));A=readtable(fullfile(out,'design_summary.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 1100 470]);tiledlayout(1,2,'TileSpacing','compact');co=[.15 .35 .65;.85 .32 .10];
for panel=1:2
 nexttile;hold on;
 for i=1:height(T)
  di=find(strcmp(D.DesignID,T.DesignID{i}));color=co(1+strcmp(T.Type{i},'PIDF'),:);mark='o';if strcmp(T.Type{i},'PIDF'),mark='^';end
  if panel==1,y=T.Recovery_s(i);else,ai=find(strcmp(A.DesignID,T.DesignID{i}));y=A.AllTargetsPassCases(ai);end
  scatter(T.CommandNoiseGain100(i),y,42,color,mark,'filled','HandleVisibility','off');
  if D.IsOriginal(di),scatter(T.CommandNoiseGain100(i),y,140,'k','p','LineWidth',1.3,'HandleVisibility','off');end
 end
 plot(NaN,NaN,'o','Color',co(1,:),'MarkerFaceColor',co(1,:),'DisplayName','Positive-gain PI designs');
 plot(NaN,NaN,'^','Color',co(2,:),'MarkerFaceColor',co(2,:),'DisplayName','Positive-gain PIDF designs');
 plot(NaN,NaN,'kp','MarkerSize',11,'DisplayName','Retained original controllers');
 grid on;xlabel('|C S| at 100 rad/s (pu/Hz)');
 if panel==1,ylabel('Nominal recovery (s)');title('Finite grid: recovery and command-noise gain');else,ylabel('All-target passes out of 43 plants');title('Sensitivity to design choices');ylim([0 45]);end
 legend('Location','best','FontSize',8);
end
style_lfc_figure(f);
exportgraphics(f,fullfile(out,'design_tradeoff.pdf'),'ContentType','vector');exportgraphics(f,fullfile(out,'design_tradeoff.png'),'Resolution',180);savefig(f,fullfile(out,'design_tradeoff.fig'));close(f);
end
