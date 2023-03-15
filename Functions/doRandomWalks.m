function [etas, resampledDemandData, runData] = doRandomWalks(dd, runLength, loadData, binEdges, binWidth, N, keepAllRunData)
% Do N random walks and return the overlap factor
    etas = [];
    runData = {};
    for n = 1:N
        tic;
        dd.randSeed = []; % Make sure the rng seed is cleared
        % Do a random walk of the transition matrix to generate a new dataset
        stepSize = 5;
        resampledDemandData = dd.genSampleData(runLength*3/stepSize, dd.rawData.datetime_utc_measured(1), 'initialMarkovStates', dd.dataBin(1,:));
        resampledDemandData = timetable(resampledDemandData(:,1),'TimeStep',minutes(stepSize),'VariableNames',{'PowerKw'});
        resampledDemandData = retime(resampledDemandData,'regular','linear','TimeStep',minutes(1));
        if keepAllRunData
            runData{end+1} = resampledDemandData;
        end
        eta = calculateEtaOverlap(loadData, binEdges, resampledDemandData, binWidth);
        etas(end+1) = eta;
        
        fprintf(1, 'Run %d - %0.2fsec\n', n, toc);
    end
end