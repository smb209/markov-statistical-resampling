function loadData = runTwoLoadCase(runId1, runId2, outdir)
    csvFile1 = sprintf('./LoadData/%d.csv',runId1);
    csvFile2 = sprintf('./LoadData/%d.csv',runId2);
    loadDataA = import5MinLoadData(char(csvFile1));
    loadDataB = import5MinLoadData(char(csvFile2));
    
    loadDataTotal = loadDataA;
    loadDataTotal.total_demand_kw = loadDataA.total_demand_kw + loadDataB.total_demand_kw;
    loadDataTotal.averaged = loadDataA.averaged + loadDataB.averaged;
%     loadDataTotal = retime(loadDataTotal,'regular','linear','TimeStep',minutes(5));
    
    fig_h = MarkovModel.plotDataWithDistribution(...
        loadDataTotal.total_demand_kw,...
        loadDataTotal.Properties.RowTimes,...
        'title',sprintf('Source Demand Data [mu: %0.1fkw]',mean(loadDataTotal.total_demand_kw)),...
        'dataName','Power Demand [kw]',...
        'fignum',390, ...
        'histBins',50);
    
    fprintf(1,'File: Combined %d,%d, meanLoad: %0.2f\n',runId1,runId2,mean(loadDataTotal.total_demand_kw))
    savefig(fig_h,fullfile(outdir,sprintf('%d_%d_orig_data_and_dist.fig',runId1,runId2)));
    saveas(fig_h,fullfile(outdir,sprintf('%d_%d_orig_data_and_dist.png',runId1,runId2)));
    loadData = loadDataTotal;
    
    loadData.net_demand_kw = loadData.total_demand_kw;
end