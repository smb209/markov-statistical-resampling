function [dd, loadData, binEdges, binWidth] = markovProcessData(loadData, dataColumnName)
    dd = MarkovModel;
    
    loadData.minutes = hour(loadData.datetime_utc_measured) * 60 + minute(loadData.datetime_utc_measured);
    duration = loadData.datetime_utc_measured(3)-loadData.datetime_utc_measured(2); % (issue with first 2 times), updated [JC]
    
    % Import the load data as a TimeTable into the MarkovModel class
    dd.importTimeTable(loadData,{dataColumnName,'minutes'},'sampleTime', duration);
    
    % Convert the continous data into bins with size [50, 144] -> [power_kw, 5-min samples/12 hrs/day]
    dd.createStateBins('histBins',[60 142]);
%     dd.createStateBins();
    
    % Save these variables off for later
    binEdges = dd.dataBinEdges{1};
    binWidth = median(diff(dd.dataBinEdges{1}));
    
    % Flag any holes in the input data that are larger than the given step
    % size, ensure they dont contribute to the expectation matrix
    dd.computeStateValidityMap(days(10)); % JC - updated
    
    % Process the raw data to generate the Markov Transition Matrix
    dd.computeNdMarkovMatrix();
    
    % Plot the transition matrix with the interactive plot
    dd.plot2dStates('fignum',398);
    
    % Look for any terminal states and remove them, rescale the
    % probabilities of ingress nodes to normalize to 1
    dd.removeTerminalStates();
    dd.plot2dStates('fignum',399);

    % Create the expectation matrix for all the data
    dd.generateFullExpectationMatrix();
end