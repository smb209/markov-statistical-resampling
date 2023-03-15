function eta = calculateEtaOverlap(loadData, binEdges, resampledDemandData, binWidth)
    % Calculates the overlap of 2 distributions
    hd1 = histcounts(loadData.total_demand_kw,binEdges, 'normalization','pdf');
    hr1 = histcounts(resampledDemandData.PowerKw,binEdges, 'normalization','pdf');
    
    % compute overlapping index (see Pastore and Calcagni (2019))
    N = length(hd1);
    eta = 0;
    for n = 1:N
        eta = eta + min(hd1(n),hr1(n))*binWidth;
    end
end