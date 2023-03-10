% Description: this file
%  Run the full study
% clear dd; clc; close all; clear all;
warning('on');
outdir = 'outputs';

%% Run length, step size - Minutes
Ndays = 60;    % number of days in study, added [JC]
Day2Hour = 24;  % day to hour conversion, added [JC]
Hour2Min = 60;  % hour to min conversion, added [JC]
runLength = Ndays*Day2Hour*Hour2Min; % updated [JC]
dataColumnName = 'net_demand_kw';

%% Import PV Data
% load measured irradiance data
csv_file = "2017.csv";
csv_file_path = fullfile('PVData',char(csv_file));

[GHITable] = ImportGHIData(csv_file_path);

%% Case 1 - Single building load
% Load 1 Load Profiles
caseName = 'one_load_case_no_pv';
runId1 = 78;
loadData = OneLoadCase(runId1, outdir);

% Initialize the MarkovModel class and preprocess the load data into states
[dd, loadData, binEdges, binWidth] = MarkovProcessData(loadData, dataColumnName);

%% Do random walks and create plots
[etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1);
resampledDemandDataForPlotX = resampledDemandData.Time+loadData.datetime_utc_measured(1);
resampledDemandDataForPlotY = resampledDemandData.PowerKw;
dateFilter = resampledDemandDataForPlotX <= loadData.datetime_utc_measured(end);

% Plot time series for case
fignum = 901;
origX = loadData.datetime_utc_measured;
origY = loadData.total_demand_kw;
newX = resampledDemandDataForPlotX(dateFilter);
newY = resampledDemandDataForPlotY(dateFilter);

PlotRunTimeseries(fignum, origX, origY, newX, newY);

% Plot Histograms for case
binEdges = dd.dataBinEdges{1};
origY = loadData.total_demand_kw;
newY = resampledDemandData.PowerKw;
eta = etas(1);
fignum = 902;

PlotRunHistograms(fignum, origY, binEdges, newY, eta);

% Resample the data set and calculate the eta overlap for 1000 runs

% SAVED FOR SPEED % 
[etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1000);
save(sprintf('results_%s.mat',caseName),'etas')
% etas = load('results_one_load_case_no_pv.mat')

% Plot Etas
fignum = 403;
saveTikz = false;

PlotEtas(fignum, etas, saveTikz, caseName);

%% Case 2 - Two building load
% Load 2 Load Profiles
caseName = 'two_load_case_no_pv';
runId1 = 101;
runId2 = 51;
loadData = TwoLoadCase(runId1, runId2, outdir);

% Initialize the MarkovModel class and preprocess the load data into states
[dd, loadData, binEdges, binWidth] = MarkovProcessData(loadData, dataColumnName);

%% Do random walks and create plots
[etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1);
resampledDemandDataForPlotX = resampledDemandData.Time+loadData.datetime_utc_measured(1);
resampledDemandDataForPlotY = resampledDemandData.PowerKw;
dateFilter = resampledDemandDataForPlotX <= loadData.datetime_utc_measured(end);

% Plot time series for case
fignum = 901;
origX = loadData.datetime_utc_measured;
origY = loadData.total_demand_kw;
newX = resampledDemandDataForPlotX(dateFilter);
newY = resampledDemandDataForPlotY(dateFilter);

PlotRunTimeseries(fignum, origX, origY, newX, newY);

% Plot Histograms for case
binEdges = dd.dataBinEdges{1};
origY = loadData.total_demand_kw;
newY = resampledDemandData.PowerKw;
eta = etas(1);
fignum = 902;

PlotRunHistograms(fignum, origY, binEdges, newY, eta);

%% Resample the data set and calculate the eta overlap for 1000 runs

% SAVED FOR SPEED % 
[etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1000);
save(sprintf('results_%s.mat',caseName),'etas')
% etas = load('results_one_load_case_no_pv.mat')

% Plot Etas
fignum = 403;
saveTikz = false;

PlotEtas(fignum, etas, saveTikz, caseName);

%% Case 3 - Two building load with PV
% Load 2 Load Profiles with PV
caseName = 'two_load_case_with_pv';
runId1 = 101;
runId2 = 51;

loadData = TwoLoadPvCase(runId1, runId2, outdir, GHITable);

% Initialize the MarkovModel class and preprocess the load data into states
[dd, loadData, binEdges, binWidth] = MarkovProcessData(loadData, dataColumnName);


%% Do random walks and create plots
[etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1);
resampledDemandDataForPlotX = resampledDemandData.Time+loadData.datetime_utc_measured(1);
resampledDemandDataForPlotY = resampledDemandData.PowerKw;
dateFilter = resampledDemandDataForPlotX <= loadData.datetime_utc_measured(end);

% Plot time series for case
fignum = 901;
origX = loadData.datetime_utc_measured;
origY = loadData.total_demand_kw;
newX = resampledDemandDataForPlotX(dateFilter);
newY = resampledDemandDataForPlotY(dateFilter);

PlotRunTimeseries(fignum, origX, origY, newX, newY);

% Plot Histograms for case
binEdges = dd.dataBinEdges{1};
origY = loadData.total_demand_kw;
newY = resampledDemandData.PowerKw;
eta = etas(1);
fignum = 902;

PlotRunHistograms(fignum, origY, binEdges, newY, eta);

% PV Data Plot

fignum = 890;
origX = loadData.datetime_utc_measured;
origY = loadData.pv_array_output_kw;

PlotPvProfile(fignum, origX, origY);

%% Resample the data set and calculate the eta overlap for 1000 runs

% SAVED FOR SPEED % 
[etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1000);
save(sprintf('results_%s.mat',caseName),'etas')
% etas = load('results_one_load_case_no_pv.mat')

% Plot Etas
fignum = 403;
saveTikz = false;

PlotEtas(fignum, etas, saveTikz, caseName);


%% Helper functions

%% Case 1
function loadData = OneLoadCase(runId1, outdir)
    csvFile1 = sprintf('./LoadData/%d.csv',runId1);
    loadDataA = import5MinLoadData(char(csvFile1));
    
    loadDataTotal = loadDataA;
    loadDataTotal.total_demand_kw = loadDataA.total_demand_kw;
    loadDataTotal.averaged = loadDataA.averaged;
    loadDataTotal = retime(loadDataTotal,'regular','linear','TimeStep',minutes(5));
    
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


%% Case 2
function loadData = TwoLoadCase(runId1, runId2, outdir)
    csvFile1 = sprintf('./LoadData/%d.csv',runId1);
    csvFile2 = sprintf('./LoadData/%d.csv',runId2);
    loadDataA = import5MinLoadData(char(csvFile1));
    loadDataB = import5MinLoadData(char(csvFile2));
    
    loadDataTotal = loadDataA;
    loadDataTotal.total_demand_kw = loadDataA.total_demand_kw + loadDataB.total_demand_kw;
    loadDataTotal.averaged = loadDataA.averaged + loadDataB.averaged;
    loadDataTotal = retime(loadDataTotal,'regular','linear','TimeStep',minutes(5));
    
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


%% Case 3
function loadData = TwoLoadPvCase(runId1, runId2, outdir, GHITable)
    loadData = TwoLoadCase(runId1, runId2, outdir);

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

%% Imports a csv file with the GHI data and returns a timetable
function [GHITable] = ImportGHIData(csv_file_path)
    opts = detectImportOptions(csv_file_path);
    D = readtable(csv_file_path,opts); % <-- | date | time | GHI
    
    % Select the data columns
    dates = D{:,1};
    times = D{:,2};
    GHI   = D{:,3};
    
    %NOTE: there was an extra (redundant time) at 5/17/2017, 7:32 AM.
    indices = find(GHI(:)<0);
    GHI(indices)=0;
    
    % There are some spurious mesurements, saturate this at an average negative
    % value
    GHI(GHI < -15) = -15;
    
    % Convert PV data to a timetable
    tmpTimes = datetime(times,'InputFormat','HH:mm');
    dateTimeVec = dates + hours(hour(tmpTimes)) + minutes(minute(tmpTimes)) + years(2000);
    GHITable = timetable(dateTimeVec,GHI);
end

%% Calculates the overlap of 2 distributions
function eta = CalculateOverlap(loadData, binEdges, resampledDemandData, binWidth)
    hd1 = histcounts(loadData.total_demand_kw,binEdges, 'normalization','pdf');
    hr1 = histcounts(resampledDemandData.PowerKw,binEdges, 'normalization','pdf');
    
    % compute overlapping index (see Pastore and Calcagni (2019))
    N = length(hd1);
    eta = 0;
    for n = 1:N
        eta = eta + min(hd1(n),hr1(n))*binWidth;
    end
end

%% Do N random walks and return the overlap factor
function [etas, resampledDemandData] = DoRandomWalks(dd, runLength, loadData, binEdges, binWidth, N)
    etas = [];
    dd.randSeed = []; % Make sure the rng seed is cleared
    for n = 1:N
        
        tic;
        % Do a random walk of the transition matrix to generate a new dataset
        stepSize = 5;
        resampledDemandData = dd.genSampleData(runLength*3/stepSize);
        resampledDemandData = timetable(resampledDemandData(:,1),'TimeStep',minutes(stepSize),'VariableNames',{'PowerKw'});
        resampledDemandData = retime(resampledDemandData,'regular','linear','TimeStep',minutes(1));
        
        eta = CalculateOverlap(loadData, binEdges, resampledDemandData, binWidth);
        etas(end+1) = eta;
        
        fprintf(1, 'Run %d - %0.2fsec\n', n, toc);
    end
end

%% 
function [dd, loadData, binEdges, binWidth] = MarkovProcessData(loadData, dataColumnName)
    dd = MarkovModel;
    
    loadData.minutes = hour(loadData.datetime_utc_measured) * 60 + minute(loadData.datetime_utc_measured);
    duration = loadData.datetime_utc_measured(3)-loadData.datetime_utc_measured(2); % (issue with first 2 times), updated [JC]
    
    % Import the load data as a TimeTable into the MarkovModel class
    dd.importTimeTable(loadData,{dataColumnName,'minutes'},'sampleTime', duration);
    
    % Convert the continous data into bins with size [50, 144] -> [power_kw, 5-min samples/12 hrs/day]
    dd.createStateBins('histBins',[60 144]);
    
    % Save these variables off for later
    binEdges = dd.dataBinEdges{1};
    binWidth = median(diff(dd.dataBinEdges{1}));
    
    % Flag any holes in the input data that are larger than the given step
    % size, ensure they dont contribute to the expectation matrix
    dd.computeStateValidityMap(minutes(60)); % JC - updated
    
    % Process the raw data to generate the Markov Transition Matrix
    dd.computeNdMarkovMatrix();
    
    % Plot the transition matrix with the interactive plot
    dd.plot2dStates('fignum',398);
    
    % Look for any terminal states and remove them, rescale the
    % probabilities of ingress nodes to normalize to 1
    dd.removeTerminalStates();
    
    % Create the expectation matrix for all the data
    dd.generateFullExpectationMatrix();
end

%% Plot histograms for a given case
function fig_h = PlotRunHistograms(fignum, origY, binEdges, newY, eta)
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

%% Plot the timeseries of the original and resampled data
function fig_h = PlotRunTimeseries(fignum, origX, origY, newX, newY)
    fig_h = figure(fignum);
    clf
    ax_h = subplot(1,1,1);
    hold(ax_h,'on');
    p1 = plot(origX,origY,'DisplayName','Measured Load');
    p2 = plot(newX, newY,'DisplayName','Resampled Load');
    legend(ax_h,'show');
    set(fig_h,'Color',[1 1 1]);
    set(ax_h,'FontName','Latin Modern Math','FontSize',22,'LineWidth',2.0);
    ylabel('$Load[kw]$','FontSize',18,'Interpreter','latex')
    xlabel('$Date$','FontSize',18,'Interpreter','latex')
    box off
end

%% Plot the distribution of etas over the N runs
function fig_h = PlotEtas(fignum, etas, saveTikz, runName)
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

%% Plot the PV profile
function PlotPvProfile(fignum, origX, origY)
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