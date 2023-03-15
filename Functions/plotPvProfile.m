function plotPvProfile(fignum, origX, origY)
% Plot the PV profile
    fig_h = figure(fignum);
    clf
    ax_h = subplot(1,1,1);
    hold(ax_h,'on');
    p1 = plot(origX,origY,'DisplayName','PV Output [kw]');
    legend(ax_h,'show');
    set(fig_h,'Color',[1 1 1]);
    set(ax_h,'FontName','Latin Modern Math','FontSize',22,'LineWidth',2.0);
    ylabel('$Load[kw]$','FontSize',18,'Interpreter','latex')
    xlabel('$Date$','FontSize',18,'Interpreter','latex')
    box off
end