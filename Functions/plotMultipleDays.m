function fig_h = plotMultipleDays(fignum, runData, loadData)
fig_h = figure(fignum);
clf
ax_h = subplot(1,1,1);
hold(ax_h,'on');

for ix = 1:numel(runData)
    resampledDemandDataForPlotX = runData{ix}.Time+loadData.datetime_utc_measured(1);
    resampledDemandDataForPlotY = runData{ix}.PowerKw;
    newX = resampledDemandDataForPlotX; %(dateFilter);
    newY = resampledDemandDataForPlotY; %(dateFilter);
    plot(newX, newY,'DisplayName','Resampled Load');
end

% legend(ax_h,'show');
set(fig_h,'Color',[1 1 1]);
set(ax_h,'FontName','Latin Modern Math','FontSize',22,'LineWidth',2.0);
ylabel('$Load[kw]$','FontSize',18,'Interpreter','latex')
xlabel('$Date$','FontSize',18,'Interpreter','latex')
box off
end