function plot_nyquist_lfc(loops,controllers,out)
% Plot positive and negative frequency branches, with direction arrows.
% The origin indentation used for the winding count is outside this finite view.
w=logspace(-3,3,3000)';col=[.04 .36 .65;.85 .32 .08];
f=figure('Visible','off','Color','w','Position',[100 100 1100 520]);tiledlayout(1,2,'TileSpacing','compact');
for panel=1:2
 nexttile;hold on;
 for ci=1:2
  h=squeeze(freqresp(loops{ci},w));
  plot(real(h),imag(h),'Color',col(ci,:),'LineWidth',1.5,'DisplayName',controllers(ci).name);
  plot(real(h),-imag(h),'--','Color',col(ci,:),'LineWidth',1.1,'HandleVisibility','off');
  for wa=[1.7 3.5 7]
   [~,k]=min(abs(w-wa));v=h(k+12)-h(k);v=.14*v/abs(v);
   quiver(real(h(k)),imag(h(k)),real(v),imag(v),0,'Color',col(ci,:),'LineWidth',1.5,'MaxHeadSize',1.5,'HandleVisibility','off');
  end
 end
 plot(-1,0,'rx','MarkerSize',10,'LineWidth',2,'DisplayName','Critical point -1');
 xline(0,':','HandleVisibility','off');yline(0,':','HandleVisibility','off');grid on;axis equal;
 if panel==1,xlim([-1.8 .8]);ylim([-3 3]);title('Nyquist branches');else,xlim([-1.3 .3]);ylim([-1 1]);title('Detail near -1');end
 xlabel('Real L(j\omega)');ylabel('Imaginary L(j\omega)');legend('Location','southwest');
end
sgtitle('Solid: positive frequency; dashed: negative; arrows: increasing positive frequency','FontSize',12);
style_lfc_figure(f);
exportgraphics(f,fullfile(out,'nyquist.png'),'Resolution',180);exportgraphics(f,fullfile(out,'nyquist.pdf'),'ContentType','vector');savefig(f,fullfile(out,'nyquist.fig'));close(f);
end
