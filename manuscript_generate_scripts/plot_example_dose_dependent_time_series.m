setup_environment

%% additional options

savePlot = 1;
plotIt = 1;
chanInt = 14;
smoothPP = 1;
rerefMode = 'median';
idx = 7;


labelTotal = [];
keepsTotal = [];
awinsTotal = [];

%%
sid = SIDS{idx};

switch(sid)
    
    case 'd5cd55'
        stims = [54 62];
        rerefChans = [2:40];
        goods = sort([44 45 46 52 53 55 60 61 63]);
        betaChan = 53;
        bads = [1 49 58 59];
        
        % have to set t_min and t_max for each subject
        %t_min = 0.004833;
        % t_min of 0.005 would find the really early ones
        t_min = 0.006;
        t_max = 0.06;
        
    case 'c91479'
        stims = [55 56];
        betaChan = 64;
        goods = sort([ 39 40 47 48 63 64]);
        bads = [1 2 3 31 57];
        rerefChans = [4:29 32:37 41:45 49:52 57:61 ];
        t_min = 0.005;
        t_max = 0.036;
    case '7dbdec'
        % rerefChans = [1:3 7 15 1:16 17:19 22:24 33:56 58:64]; how it
        % was for paper
        
        rerefChans = [1:16 17:19 22:24 33:56 58:64]; % without doubling up
        stims = [11 12];
        chans = [4 5 14];
        goods = sort([4 5 10 13]);
        betaChan = 4;
        t_min = 0.007;
        t_max = 0.048;
        bads = [8 57];
        
    case '9ab7ab'
        stims = [59 60];
        betaChan = 51;
        rerefChans = [2:8 10:40 45:48 54:56];
        goods = sort([42 43 49 50 51 52 53 57 58]);
        t_min = 0.006;
        t_max = 0.06;
        bads = [1 9 10 35 43];
        
    case '702d24'
        rerefChans = [1:4 6:12 15:22 24 25:27 33:40 41:43 45:51 53:58 62:64];
        betaChan = 5;
        stims = [13 14];
        goods = [ 5 ];
        t_min = 0.008;
        t_max = 0.046;
        bads = [23 27 28 29 30 32 44 52 60];
        
    case 'ecb43e' % added DJC 7-23-2015
        rerefChans = [1:40 41:44 49:52];
        stims = [56 64];
        betaChan = 55;
        goods = sort([55 63 54 47 48]);
        bads = [57:64];
        
        t_min = 0.006;
        t_max = 0.06;
    case '0b5a2e' % added DJC 7-23-2015
        rerefChans = [1:8 9:12 17:20 24 25:28 33:37 38 41:48 49:64];
        stims = [22 30];
        betaChan = 23;
        goods = [14 21 23 31];
        bads = [24 25 29];
        t_min = 0.005;
        t_max = 0.06;
    case '0b5a2ePlayback' % added DJC 7-23-2015
        rerefChans = [1:8 9:12 17:20 24 25:28 33:37 38 41:48 49:64];
        stims = [22 30];
        betaChan = 31;
        goods = sort([12 13 14 15 16 21 23 31 32 39 40]);
        bads = [24 25 29];
        t_min = 0.005;
        t_max = 0.06;
    otherwise
        error('unknown SID entered');
end

badsTotal = [stims bads];

chans = [1:64];
chans(ismember(chans, badsTotal)) = [];
%% load in the trigger data

if strcmp(sid,'0b5a2ePlayback')
    load(fullfile(folderTiming, ['0b5a2e_tables.mat']), 'bursts', 'fs', 'stims');
    delay = 577869;
else
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
delayDelivery = 14;

% add 14 samples to approximately account for the delay between the
% stimulation trigger command and actual registration
stims(2,:) = stims(2,:) + delayDelivery;
bursts(2,:) = bursts(2,:) + delayDelivery;
bursts(3,:) = bursts(3,:) + delayDelivery;

%figure
%% process each ecog channel individually

sigChans = {};
shuffleChans = {};
CCEPbyNumStim = {};

dataForAnova = {};

ZscoredDataForAnova = {};

%% do referencing on list of channels

for chan = rerefChans
    
    %% load in ecog data for that channel
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
    else
        error 'unknown sid';
    end
    
    presamps = round(0.05*efs);
    postsamps = round(0.120*efs);
    
    ptis = round(stims(2,pts)/fac);
    
    t = (-presamps:postsamps)/efs;
    
    
    wins = squeeze(getEpochSignal(eco', ptis-presamps, ptis+postsamps+1));
    winsReref(:,:,chan) = wins;
end

switch(rerefMode)
    case 'mean'
        rerefQuant = mean(winsReref(:,:,rerefChans),3);
        
    case 'median'
        rerefQuant = median(winsReref(:,:,rerefChans),3);
end

%% now do peak to peak
% set statistical threshold

for chan = chanInt
    fprintf('loading in ecog data for %s:\n',sid);
    fprintf('channel %d:\n',chan);
    tic;
    
    grp = floor((chan-1)/16);
    ev = sprintf('ECO%d',grp+1);
    achan = chan - grp*16;
    
    load(fullfile(folderECoGData,[sid '_ECoG.mat']),ev);
    dataStruct = eval(ev);
    
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
        
    else
        error 'unknown sid';
    end
    
    presamps = round(0.05*efs);
    postsamps = round(0.120*efs);
    
    ptis = round(stims(2,pts)/fac);
    
    t = (-presamps:postsamps)/efs;
    
    wins = squeeze(getEpochSignal(eco', ptis-presamps, ptis+postsamps+1));
    wins = wins - rerefQuant;
    pstims = stims(:,pts);
    
    % considered a baseline if it's been at least N seconds since the last
    % burst ended
    
    baselines = pstims(5,:) > 2 * fs;
    
    if (sum(baselines) < 100)
        warning('N baselines = %d.', sum(baselines));
    end
    
    types = unique(bursts(5,pstims(4,:)));
    
    %DJC - modify suffix to list conditioning type
    suffix = arrayfun(@(x) num2str(x), types, 'uniformoutput', false);
    %
    suffix = cell(1,4);
    suffix{1} = 'Phase 1';
    suffix{2} = 'Phase 2';
    suffix{3} = 'null condition';
    suffix{4} = 'Phase 3';
    
    nullType = 2;
    maxIndex = min(length(types),2);
    
    for typei = 1:maxIndex
        awins = wins-repmat(mean(wins(t<-0.005 & t>-0.05,:),1), [size(wins, 1), 1]);
        
        probes = pstims(5,:) < .5*fs & bursts(5,pstims(4,:))==types(typei);
        
        if (sum(probes) < 100)
            warning('N probes = %d.', sum(probes));
        end
        
        
        if (types(typei) == nullType)
            label = nan(1,length(bursts(4,pstims(4,:))));
            label(baselines) = 0;
            label(probes) = 1 ;
            
        elseif (types(typei) ~= nullType)
            label = bursts(4,pstims(4,:));
            label(baselines) = 0;
            labelGroupStarts = [1 3 5];
            labelGroupEnds   = [labelGroupStarts(2:end) Inf];
            
            for gIdx = 1:length(labelGroupStarts)
                labeli = label >= labelGroupStarts(gIdx) & label < labelGroupEnds(gIdx);
                label(labeli) = gIdx;
            end
            
        end
        
        keeps = probes | baselines;
        load('line_colormap.mat');
        
        kwins = awins(:, keeps);
        klabel = label(keeps);
        ulabels = unique(klabel);
        colors = cm(round(linspace(1, size(cm, 1), length(ulabels))), :);
        
        tBegin = t_min;
        tEnd = t_max;
        
        labelTotal = [labelTotal label];
        keepsTotal = [keepsTotal keeps];
        awinsTotal = [awinsTotal awins];
    end
    
    keepsTotal = logical(keepsTotal);
    %%
    figure
    % this sets the figure to be the whole screen
    set(gcf, 'Units', 'inches', 'Position', [0 0 8 4]);
    prettylineNoStd(1e3*t,1e6*awinsTotal(:, keepsTotal), labelTotal(keepsTotal), colors);
    xlim(1e3*[min(t) max(t)]);
    yl = ylim;
    yl(1) = min(-10, max(yl(1),-140*4));
    yl(2) = max(10, min(yl(2),100*4));
    
    yl(1) = min(-10, max(yl(1),-340*4));
    yl(2) = max(10, min(yl(2),300*4));
    
    ylim(yl);
    
    xlim([-2.5 60])
    ylim([-300 300])
    highlight(gca, [0 t_min*1e3], [], [.5 .5 .5]) %this is the part that plots that stim window
    
    %  vline(1e3*7/efs);
    vline(0);
    obj = scalebar;
    obj.XLen = 30;              %X-Length, 10.
    obj.XUnit = 'ms';            %X-Unit, 'm'.
    obj.YLen = 200;
    obj.YUnit = '\muV';
    
    obj.Position = [20,-130];
    obj.hTextX_Pos = [5,-30]; %move only the LABEL position
    obj.hTextY_Pos =  [30,-15];
    obj.hLineY(2).LineWidth = 5;
    obj.hLineY(1).LineWidth = 5;
    obj.hLineX(2).LineWidth = 5;
    obj.hLineX(1).LineWidth = 5;
    obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'    if savePlot
    
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    %set(gca,'visible','off')
    
    for index = 0:max(labelTotal) 
        awinsSelect = awinsTotal(:,(labelTotal == index)& keepsTotal);
        awinsSelectMean = mean (awinsSelect,2);
        [signalPP,pkLocs,trLocs] =  extract_PP_betaStim(awinsSelectMean,t,t_min,t_max,smoothPP);
    end
    
    %     SaveFig(folderPlots, sprintf(['EP-phase-%d-sid-%s-chan-%d'],typei,sid, chan,type,signalType), 'svg');
    SaveFig(folderPlots, sprintf(['EP-example-time-series-sid-%s-chan-%d'],sid, chan), 'png','-r600');
    SaveFig(folderPlots, sprintf(['EP-example-time-series-sid-%s-chan-%d'],sid, chan), 'eps','-r600');
    
    
    % end
end







