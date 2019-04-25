locationsDir = fileparts(which('master_script_betaStim'));
SIDS = {'d5cd55', 'c91479', '7dbdec', '9ab7ab', '702d24', 'ecb43e','0b5a2e','0b5a2ePlayback'};
prefixDirectory = locationsDir;

folderTiming = fullfile(prefixDirectory,'data','stim_timing_data');
folderECoGData = fullfile(prefixDirectory,'data','ECoG_data');
folderCoords = fullfile(prefixDirectory,'coordinates');
folderPlots = fullfile(prefixDirectory,'output_plots');
folderTestFilter = fullfile(prefixDirectory,'test_real_time_filter');