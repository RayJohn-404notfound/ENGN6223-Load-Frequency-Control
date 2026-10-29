function style_lfc_figure(f)
%STYLE_LFC_FIGURE Explicit print styling, independent of desktop theme.
set(f,'Color','w');
ax=findall(f,'Type','axes');set(ax,'Color','w','XColor','k','YColor','k','GridColor',[.65 .65 .65]);
set(findall(f,'Type','text'),'Color','k');
set(findall(f,'Type','legend'),'Color','w','TextColor','k','EdgeColor',[.6 .6 .6]);
end
