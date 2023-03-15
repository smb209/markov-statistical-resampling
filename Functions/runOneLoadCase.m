function loadData = runOneLoadCase(runId1, outdir)
    csvFile1 = sprintf('./LoadData/%d.csv',runId1);
    loadDataA = import5MinLoadData(char(csvFile1));
    
    loadDataTotal = loadDataA;
    loadDataTotal.total_demand_kw = loadDataA.total_demand_kw;
    loadDataTotal.averaged = loadDataA.averaged;
%     loadDataTotal = retime(loadDataTotal,'regular','nearest','TimeStep',minutes(5));
    
    fig_h = MarkovModel.plotDataWithDistribution(...
        loadDataTotal.total_demand_kw,...
        loadDataTotal.Properties.RowTimes,...
        'title',sprintf('Source Demand Data [mu: %0.1fkw]',mean(loadDataTotal.total_demand_kw)),...
        'dataName','Power Demand [kw]',...
        'fignum',390, ...
        'histBins',50);
    
    fprintf(1,'File: Run %d meanLoad: %0.2f\n',runId1,mean(loadDataTotal.total_demand_kw))
    savefig(fig_h,fullfile(outdir,sprintf('%d_orig_data_and_dist.fig',runId1)));
    saveas(fig_h,fullfile(outdir,sprintf('%d_orig_data_and_dist.png',runId1)));
    loadData = loadDataTotal;
    
    loadData.net_demand_kw = loadData.total_demand_kw;
end