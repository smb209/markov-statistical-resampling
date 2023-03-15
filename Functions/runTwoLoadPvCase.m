function loadData = runTwoLoadPvCase(runId1, runId2, outdir, GHITable)
    loadData = runTwoLoadCase(runId1, runId2, outdir);

    % Correlate the PV Data with the load data
    pvYear = year(GHITable.dateTimeVec(1));
    timeStart = min(loadData.datetime_utc_measured);
    timeEnd = max(loadData.datetime_utc_measured);
    
    % Create a timerange from the month and day to select the appropriate GHI
    dtStart = datetime(pvYear,month(timeStart),day(timeStart));
    dtEnd = datetime(pvYear,month(timeEnd),day(timeEnd));
    tr = timerange(dtStart,dtEnd);
    
    GHIForLoad = GHITable(tr,:);
    
    % Align the GHI samples with the load samples so we can sum the vectors
    dateDiff = loadData.datetime_utc_measured(1) - GHIForLoad.dateTimeVec(1);
    GHIForLoad2 = GHIForLoad;
    GHIForLoad2.dateTimeVec = GHIForLoad2.dateTimeVec + dateDiff;
    GHIForLoad2 = retime(GHIForLoad2,loadData.datetime_utc_measured,"nearest");
    
    % Select a size of the PV array and scale the GHI data to that size using the max recorded GHI
    maxGHI = max(GHITable.GHI);
    
    % Array size is 10%+ of the max load
    arraySizeKW = ceil(max(loadData.total_demand_kw)*0.1);
    scaleFactor = arraySizeKW / maxGHI;
    
    % Assign to both timetables for now
    GHIForLoad2.array_output_kw = GHIForLoad2.GHI * scaleFactor;
    loadData.pv_array_output_kw = GHIForLoad2.GHI * scaleFactor;
    loadData.net_demand_kw = loadData.total_demand_kw - loadData.pv_array_output_kw;
end