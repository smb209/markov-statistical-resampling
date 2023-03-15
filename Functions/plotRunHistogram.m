function fig_h = plotRunHistogram(fignum, origY, binEdges, newY, eta)
% Plot histograms for a given case
    fig_h = figure(fignum);
    clf
    hold on
    ax2 = subplot(1,1,1);
    hd = histogram(ax2,origY,binEdges, 'normalization','pdf');
    hr = histogram(ax2,newY,binEdges, 'normalization','pdf');
    
    set(hd,'FaceColor','black','EdgeColor',[.5 .5 .5],'FaceAlpha',.3)
    set(hr,'FaceColor','black','EdgeColor',[.5 .5 .5],'FaceAlpha',.3)
    ylabel(ax2,'pdf$(x)$, pdf$(\hat{x})$', 'interpreter','latex')
    xlabel(ax2,'$x\,$[W]', 'interpreter','latex')
    set(ax2,'Fontsize', 20,'Fontname','Times');
    set(ax2,'color','white');
    
    % Dont include titles in the matlab plot and add it later in tex
    title(sprintf('Empirical vs. Resampled Probability Density Functions\n$eta=%0.1f\\%%$',eta*100),'interpreter','latex')
    % %set(gca, 'XDir','reverse')
    text(0,1,sprintf('$eta=%0.1f\\%%$',eta*100),'interpreter','latex')
    set(gcf,'Color',[1 1 1]);
end