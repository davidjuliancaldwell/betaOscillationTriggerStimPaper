locationsDir = pwd;
setenv('subject_dir',fullfile(locationsDir,'coordinates'));
folderData = fullfile(locationsDir,'data');
folderECoGData = fullfile(folderData,'ECoG_data');
folderTiming = fullfile(folderData,'stim_timing_data');
folderPhase = fullfile(folderData,'phase_data');
folderEP = fullfile(folderData,'EP_data');
folderCoords = fullfile(locationsDir,'coordinates');
folderPlots = fullfile(locationsDir,'plots');
folderTestFilter = fullfile(locationsDir,'test_real_time_filter');