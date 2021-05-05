figure;
ax(1) = axes('Units','normalized',...
    'Position',[.075, .075, .3, .35]);
ax(2) = axes('Units','normalized',...
    'Position',[.525, .075, .3, .35]);
ax(3) = axes('Units','normalized',...
    'Position',[.075, .575, .3, .35]);
ax(4) = axes('Units','normalized',...
    'Position',[.525, .575, .3, .35]);
for n =1:4
    imagesc(10^n * rand(10),'Parent',ax(n))
    set(gcf,'Position', [0 0 1000 1000]);
    set(gca,'FontSize', 50, 'LineWidth', 1);
    pbaspect([1 1 1]);
    c = colorbar('peer',ax(n));
    set(c,'Units','normalized');
    pC = get(c,'Position');
    pAx = get(ax(n),'Position');
    set(c,'Position',[pAx(1)+pAx(3)+.075, pC(2), .05, pC(4)])
end