%% script to look at burst hists
% written by DJC 1-8-2015

setup_environment
%sid = input('what is is the subject sid\n','s');
for sidInd = 1:1
    sid = SIDS{sidInd};
    switch(sid)
        case 'd5cd55'
            subject_num = '1';
            typeCell = {'180'};
            
        case 'c91479'
            subject_num = '2';
            typeCell = {'180','0'};
            
        case '7dbdec'
            subject_num = '3';
            typeCell = {'180'};
            
            
        case '9ab7ab'
            subject_num = '4';
            typeCell = {'270'};
            
        case '702d24'
            subject_num = '5';
            typeCell = {'270','90'};
            
        case 'ecb43e' % added DJC 7-23-2015
            subject_num = '6';
            typeCell = {'270','90','Null','Random'};
            
        case '0b5a2e' % added DJC 7-23-2015
            subject_num = '7';
            typeCell = {'270','90','Null'};
            
        case '0b5a2ePlayback' % added DJC 7-23-2015
            
            typeCell = {'270','90','Null'};
    end
    load([sid '_stim_table_data_file.mat'])
    
    %%
    
    if strcmp(sid,'0b5a2ePlayback')
        load(fullfile(folderTiming, ['0b5a2e_tables.mat']), 'bursts', 'fs', 'stims');
    else % for other subjects ( not the last 2)
        load(fullfile(folderTiming, [sid '_tables.mat']), 'bursts', 'fs', 'stims');
    end
    % drop any stims that happen in the first 500 milliseconds
    stims(:,stims(2,:) < fs/2) = [];
    
    % drop any probe stimuli without a corresponding pre-burst/post-burst
    bads = stims(3,:) == 0 & (isnan(stims(4,:)) | isnan(stims(6,:)));
    stims(:, bads) = [];
    
    % adjust stim and burst tables for 0b5a2e playback case
    
    if strcmp(sid,'0b5a2ePlayback')
        stims(2,:) = stims(2,:)+delay;
        bursts(2,:) = bursts(2,:) + delay;
        bursts(3,:) = bursts(3,:) + delay;
    end
    
    % get rid of d5cd55 bursts at beginning?
    if strcmp(sid,'1')
        bursts = bursts(:,(bursts(3,:)>36536266));
    end
    
    %%
    burst_hist(subject_num,bursts,typeCell,folderPlots)
    %burst_timing(sid,bursts)
    
end
