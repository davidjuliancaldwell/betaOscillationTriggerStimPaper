locationsDir = fileparts(which('master_script_betaStim'));
SIDS = {'d5cd55', 'c91479', '7dbdec', '9ab7ab', '702d24', 'ecb43e','0b5a2e','0b5a2ePlayback'};
setenv('subject_dir',fullfile(locationsDir,'coordinates'));

folderData = fullfile(locationsDir,'data');
folderECoGData = fullfile(folderData,'ECoG_data');
folderTiming = fullfile(folderData,'stim_timing_data');
folderPhase = fullfile(folderData,'phase_data');
folderEP = fullfile(folderData,'EP_data');
folderCoords = fullfile(locationsDir,'coordinates');
folderPlots = fullfile(locationsDir,'output_plots');
folderTestFilter = fullfile(locationsDir,'test_real_time_filter');
folderOutput = fullfile(locationsDir,'output_table');

