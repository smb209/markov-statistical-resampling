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

[GHITable] = importGhiData(csv_file_path);

%% Case 1 - Single building load
% Load 1 Load Profiles
caseName = 'one_load_case_no_pv';
runId1 = 78;
loadData = runOneLoadCase(runId1, outdir);

doRunAnalysis(loadData, dataColumnName, runLength, caseName, true);

%% Case 2 - Two building load
% Load 2 Load Profiles
caseName = 'two_load_case_no_pv';
runId1 = 101;
runId2 = 51;
loadData = runTwoLoadCase(runId1, runId2, outdir);

doRunAnalysis(loadData, dataColumnName, runLength, caseName, true);


%% Case 3 - Two building load with PV
% Load 2 Load Profiles with PV
caseName = 'two_load_case_with_pv';
runId1 = 101;
runId2 = 51;

loadData = runTwoLoadPvCase(runId1, runId2, outdir, GHITable);

doRunAnalysis(loadData, dataColumnName, runLength, caseName, true);

%% PV Data Plot

fignum = 890;
origX = loadData.datetime_utc_measured;
origY = loadData.pv_array_output_kw;

plotPvProfile(fignum, origX, origY);


%% Helper Function
function doRunAnalysis(loadData, dataColumnName, runLength, caseName, runAll)
% Initialize the MarkovModel class and preprocess the load data into states
[dd, loadData, binEdges, binWidth] = markovProcessData(loadData, dataColumnName);


%% Do random walks and create plots
[etas, resampledDemandData] = doRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1, false);
resampledDemandDataForPlotX = resampledDemandData.Time+loadData.datetime_utc_measured(1);
resampledDemandDataForPlotY = resampledDemandData.PowerKw;
dateFilter = resampledDemandDataForPlotX <= loadData.datetime_utc_measured(end);

% Plot time series for case
fignum = 901;
origX = loadData.datetime_utc_measured;
origY = loadData.total_demand_kw;
newX = resampledDemandDataForPlotX(dateFilter);
newY = resampledDemandDataForPlotY(dateFilter);

fig_h = plotRunTimeseries(fignum, origX, origY, newX, newY);
savefig(fig_h,sprintf('plt_%s_%s.fig',caseName,'runTimeSeries'));

% Plot Histograms for case
binEdges = dd.dataBinEdges{1};
origY = loadData.total_demand_kw;
newY = resampledDemandData.PowerKw;
eta = etas(1);
fignum = 902;

fig_h = plotRunHistogram(fignum, origY, binEdges, newY, eta);
savefig(fig_h,sprintf('plt_%s_%s.fig',caseName,'oneRunHistogram'));

%% Resample the data set and calculate the eta overlap for 1000 runs
%
% Bypassed % 
if (runAll)
    [etas, resampledDemandData] = doRandomWalks(dd, runLength, loadData, binEdges, binWidth, 1000, false);
    save(sprintf('results_%s.mat',caseName),'etas')
    % etas = load('results_one_load_case_no_pv.mat')
    
    % Plot Etas
    fignum = 403;
    saveTikz = false;
    
    fig_h = plotEtaDistribution(fignum, etas, saveTikz, caseName);
    savefig(fig_h,sprintf('plt_%s_%s.fig',caseName,'etasN1000'));
end

%%  Do 10 resamples
fignum = 909;

singleDayRunLength = 1*14*60/3.5-10;
[~, ~, runData] = doRandomWalks(dd, singleDayRunLength, loadData, binEdges, binWidth, 50, true);
fig_h = plotMultipleDays(fignum, runData, loadData);
savefig(fig_h,sprintf('plt_%s_%s.fig',caseName,'superimposedSynthetic'));

%%
fig_h = plotSuperimposedLoadData(910,loadData.datetime_utc_measured, loadData.net_demand_kw);
savefig(fig_h,sprintf('plt_%s_%s.fig',caseName,'superimposedLoadData'));
end