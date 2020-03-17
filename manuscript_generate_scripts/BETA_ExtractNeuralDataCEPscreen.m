%% Extract neural data for screening CEPs.

% Constants
%close all;clear all;clc
saveIt = 0;
SUB_DIR = fullfile(myGetenv('subject_dir'));
cd(fileparts(which('BETA_ExtractNeuralDataCEPscreen')));
locationsDir = pwd;
META_DIR = fullfile(locationsDir, '..','data','stim_timing_data');
%folderCoords = fullfile(locationsDir,'..','coordinates');
%%
for idx = 1:7
    sid = SIDS{idx};
    
    switch(sid)
        
        case 'd5cd55'
            stimChans = [54 62];
            betaChan = 53;
            rerefChans = [2:40];
            
        case 'c91479'
            stims = [55 56];
            stimChans = [55 56];
            betaChan = 64;
            rerefChans = [4:29 32:37 41:45 49:52 57:61 ];
            
        case '7dbdec'
            stimChans = [11 12];
            betaChan = 4;
            rerefChans = [1:16 17:19 22:24 33:56 58:64]; % without doubling up
            
        case '9ab7ab'
            stimChans = [59 60];
            betaChan = 51;
            rerefChans = [2:8 10:40 45:48 54:56];
            
        case '702d24'
            stimChans = [13 14];
            betaChan = 5;
            rerefChans = [1:4 6:12 15:22 24 25:27 33:40 41:43 45:51 53:58 62:64];
            
        case 'ecb43e'
            block = 'BetaPhase-3';
            stimChans = [56 64];
            betaChan = 55;
            rerefChans = [1:40 41:44 49:52];
            
            %         chans = [47 55]; want to look at all channels
        case '0b5a2e' % added DJC 7-23-2015
            stimChans = [22 30];       
            betaChan = 31;
            rerefChans = [1:8 9:12 17:20 24 25:28 33:37 38 41:48 49:64];
            
        case '0b5a2ePlayback' % added DJC 7-23-2015
            stimChans = [22 30];
            rerefChans = [1:8 9:12 17:20 24 25:28 33:37 38 41:48 49:64];
            
    end
    chans = [1:64];
    
    %% load in the trigger data
    if strcmp(sid,'0b5a2ePlayback')
        load(fullfile(META_DIR, ['0b5a2e' '_tables.mat']), 'bursts', 'fs', 'stims');
        delay = 577869;
    else
        load(fullfile(META_DIR, [sid '_tables.mat']), 'bursts', 'fs', 'stims');
    end
    
    % drop any stims that happen in the first 500 milliseconds
    stims(:,stims(2,:) < fs/2) = [];
    
    % drop any probe stimuli without a corresponding pre-burst/post-burst
    % still want to do this for selecting conditioning - DJC, as this in the
    % next step ensures we only select conditioning up to last one before beta
    % stim train    bads = stims(3,:) == 0 & (isnan(stims(4,:)) | isnan(stims(6,:)));
    bads = stims(3,:) == 0 & (isnan(stims(4,:)) | isnan(stims(6,:)));
    stims(:, bads) = [];
    
    % adjust stim and burst tables for 0b5a2e playback case
    
    if strcmp(sid,'0b5a2ePlayback')
        
        stims(2,:) = stims(2,:)+delay;
        bursts(2,:) = bursts(2,:) + delay;
        bursts(3,:) = bursts(3,:) + delay;
        
    end
    
    % drop any stims that happen in the first 500 milliseconds
    stims(:,stims(2,:) < fs/2) = [];
    
    delayDelivery = 14;
    
    % add 14 samples to approximately account for the delay between the
    % stimulation trigger command and actual registration
    stims(2,:) = stims(2,:) + delayDelivery;
    bursts(2,:) = bursts(2,:) + delayDelivery;
    bursts(3,:) = bursts(3,:) + delayDelivery;
    
    bads = stims(3,:) == 0 & (isnan(stims(4,:)) | isnan(stims(6,:)));
    stims(:, bads) = [];
    
    %% do referencing on list of channels
    
    for chan = rerefChans
        
        % load in ecog data for that channel
        fprintf('loading in ecog data for %s:\n',sid);
        fprintf('channel %d:\n',chan);
        tic;
        
        grp = floor((chan-1)/16);
        ev = sprintf('ECO%d',grp+1);
        achan = chan - grp*16;
        
        if achan==1 || achan == 2 || achan == 4 || achan == 6
            load(fullfile(folderECoGData,[sid '_ECoG.mat']),ev);
            dataStruct = eval(ev);
        end
        eco = dataStruct.data(:,achan);
        eco = 4*eco';
        efs = dataStruct.info.SamplingRateHz;
        toc;
        
        fac = fs/efs;
        %% process triggers
        
        if (strcmp(sid, '8adc5c'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, 'd5cd55'))
            %         pts = stims(3,:)==0 & (stims(2,:) > 4.5e6);
            pts = stims(3,:)==0 & (stims(2,:) > 4.5e6) & (stims(2, :) > 36536266);
        elseif (strcmp(sid, 'c91479'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, '7dbdec'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, '9ab7ab'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, '702d24'))
            pts = stims(3,:)==0;
            %modified DJC 7-27-2015
        elseif (strcmp(sid, 'ecb43e'))
            pts = stims(3,:) == 0;
        elseif (strcmp(sid, '0b5a2e'))
            pts = stims(3,:) == 0;
        elseif (strcmp(sid, '0b5a2ePlayback'))
            pts = stims(3,:) == 0;
        elseif (strcmp(sid,'3f2113'))
            pts = stims(3,:) == 0;
            
        else
            error 'unknown sid';
        end
        
        presamps = round(0.100*efs);
        postsamps = round(0.120*efs);
        
        ptis = round(stims(2,pts)/fac);
        
        t = (-presamps:postsamps)/efs;
        
        
        wins = squeeze(getEpochSignal(eco', ptis-presamps, ptis+postsamps+1));
        winsReref(:,:,chan) = wins;
    end
    
    rerefMode = 'median';
    switch(rerefMode)
        case 'mean'
            rerefQuant = mean(winsReref(:,:,rerefChans),3);
            
        case 'median'
            rerefQuant = median(winsReref(:,:,rerefChans),3);
    end
    
    %% process each ecog channel individually
    
    for chan = chans
        %% load in ecog data for that channel
        fprintf('loading in ecog data for %s:\n',sid);
        fprintf('channel %d:\n',chan);
        tic;
        
        grp = floor((chan-1)/16);
        ev = sprintf('ECO%d',grp+1);
        achan = chan - grp*16;
        
        if achan==1 || achan == 2
            load(fullfile(folderECoGData,[sid '_ECoG.mat']),ev);
            dataStruct = eval(ev);
        end
        eco = dataStruct.data(:,achan);
        eco = 4*eco';
        efs = dataStruct.info.SamplingRateHz;
        toc;
        
        fac = fs/efs;        
        %% process triggers
        
        if (strcmp(sid, '8adc5c'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, 'd5cd55'))
            %         pts = stims(3,:)==0 & (stims(2,:) > 4.5e6);
            pts = stims(3,:)==0 & (stims(2,:) > 4.5e6) & (stims(2, :) > 36536266);
        elseif (strcmp(sid, 'c91479'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, '7dbdec'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, '9ab7ab'))
            pts = stims(3,:)==0;
        elseif (strcmp(sid, '702d24'))
            pts = stims(3,:)==0;
            %modified DJC 7-27-2015
        elseif (strcmp(sid, 'ecb43e'))
            pts = stims(3,:) == 0;
        elseif (strcmp(sid, '0b5a2e'))
            pts = stims(3,:) == 0;
        elseif (strcmp(sid, '0b5a2ePlayback'))
            pts = stims(3,:) == 0;
        elseif (strcmp(sid,'3f2113'))
            pts = stims(3,:) == 0;
            
        else
            error 'unknown sid';
        end
        
        ptis = round(stims(2,pts)/fac);
        
        % change presamps and post samps to be what Kurt wanted to look at
        presamps = round(0.1 * efs); % pre time in sec
        postsamps = round(0.12 * efs); % post time in sec
        
        t = (-presamps:postsamps)/efs;
        wins = squeeze(getEpochSignal(eco', ptis-presamps, ptis+postsamps+1));
        
        wins = wins - rerefQuant;
        awins = wins-repmat(mean(wins(t<-0.005,:),1), [size(wins, 1), 1]);
        pstims = stims(:,pts);
        
        if ~strcmp('0a80cf',sid)
            types = unique(bursts(5,pstims(4,:)));
        end
        
        % using idea of baselines
        
        baselines = pstims(5,:) > 2 * fs;
        keeper = baselines;
        
        if ~strcmp('0a80cf',sid)
            types = unique(bursts(5,pstims(4,:)));
        end
                
        kwins = awins(:,keeper);
        kwins = awins;
        clear awins
   
        ECoGData(:,:,chan) = kwins;
        clear kwins
    end
    
    ECoGDataAverage = squeeze(mean(ECoGData,2));
    
    if saveIt
        save(fullfile(OUTPUT_DIR, [sid '_baselineCCEPs.mat']), 't','ECoGData','ECoGDataAverage','-v7.3');
    end
   
    %%
    
    smallMultiples(ECoGDataAverage,t,'type1',stimChans,'type2',betaChan,'average',1);
    
    clearvars winsReref ECoGData
    
end
