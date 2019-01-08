%% script to recreate paper analyses
%Dose Dependent Enhancement of Cortically Evoked Potentials During
%Beta-Oscillation Phase Triggered Direct Cortical Stimulation of Human
%Cortex
% David.J.Caldwell

cd(fileparts(which('master_script_betaStim')));
locationsDir = pwd;
folderData = fullfile(locationsDir,'data');
folderCoords = fullfile(locationsDir,'coordinates');
folderPlots = fullfile(locationsDir,'plots');

saveIt = 1;

%% build tables with stimulation locations

A_BuildStimTablesFirst6

A_BuildStimTablesSubj7

%% extract peak to peak
B_ExtractNeuralData_PP_reref

%% plot cortex

C_PlotBrains_PP