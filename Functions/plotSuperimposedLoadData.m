function fig_h = plotSuperimposedLoadData(fignum, origX, origY)
fig_h = figure(fignum);
clf
ax_h = subplot(1,1,1);
hold(ax_h,'on');

dayNum = day(origX); %minute(origX) + hour(origX) * 60;
dayNumStart = find(diff(dayNum) ~= 0) + 1;

for nIx = 2:numel(dayNumStart)
    ixs = dayNumStart(nIx-1):dayNumStart(nIx)-1;
    x = minute(origX(ixs)) + hour(origX(ixs))*60;
    y = origY(ixs);
    plot(x, y,'DisplayName','Orig Load');
end


% legend(ax_h,'show');
set(fig_h,'Color',[1 1 1]);
set(ax_h,'FontName','Latin Modern Math','FontSize',22,'LineWidth',2.0);
ylabel('$Load[kw]$','FontSize',18,'Interpreter','latex')
xlabel('$Minute$','FontSize',18,'Interpreter','latex')
box off
end