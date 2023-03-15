function fig_h = plotEtaDistribution(fignum, etas, saveTikz, runName)
% Plot the distribution of etas over the N runs
    fig_h = figure(fignum);
    ax_h = subplot(1,1,1);
    hr = histogram(ax_h,etas,'Normalization','probability');
    
    set(fig_h,'Color',[1 1 1]);
    set(hr,'FaceColor','black','EdgeColor',[.5 .5 .5],'FaceAlpha',.3)
    set(ax_h,'FontName','Latin Modern Math','FontSize',22,'LineWidth',2.0);
    ylabel('$P$','FontSize',18,'Interpreter','latex')
    xlabel('$eta\,$','FontSize',18,'Interpreter','latex')
    box off
    
    if (saveTikz)
        cleanfigure
        matlab2tikz(sprintf('%s.tex',runName),'width','5in','height','3in',...
            'extraaxisoptions',['title style={font=\Huge},'...
            'xlabel style={font=\normalsize},'...
            'ylabel style={font=\normalsize},',...
            'ticklabel style={font=\footnotesize}']);
    end
end