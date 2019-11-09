%% script to recreate paper analyses
%Dose Dependent Enhancement of Cortically Evoked Potentials During
%Beta-Oscillation Phase Triggered Direct Cortical Stimulation of Human
%Cortex
% David.J.Caldwell

cd(fileparts(which('master_script_betaStim')));

setup_environment

saveIt = 1;
generateIntermediateData = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if generateIntermediateData
    %% build tables with stimulation locations
    
    A_BuildStimTablesFirst6
    
    A_BuildStimTablesSubj7andPlayback
    
    %% extract peak to peak
    B_ExtractNeuralData_PP_reref

    multipleSubj_GLMM_script_PP

    %% extract phase
    % this can take a long time
    B_phaseCalc_allChans_processed
    
end

%% generate table for linear mixed model

multipleSubj_GLMM_script_PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot histograms of phase delivery
burstHistogram_generation

%% plot cortex with peak to peak
C_PlotBrains_PP

%% plot phases
% phase on cortex
plotPhase_distributions_cortex

% phase for each channel, each subject
plotPhase_distributions_looping_allChans

%% plot phases and EPs together
phase_vs_peak

%% supplementary example fits
plot_example_subject_phases

%% supplementary oscilloscope test
examineExampleSigs

examinerandomSubjectDelivery

extractOscope



