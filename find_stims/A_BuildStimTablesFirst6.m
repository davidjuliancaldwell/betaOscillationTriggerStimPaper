%% Constants
%close all; clear all;clc

Z_Constants;
%% Load in the trigger data

for index = 1:6
    
    % select the subject from list
    sid = SIDS{index};
    
    load([sid '_stim_table_data_file.mat'])
    %% build a burst table with the following
    % 1 - burst id
    % 2 - burst start sample
    % 3 - burst stop sample
    % 4 - nstims in burst (do last)
    % 5 - type of conditioning (0 = falling, 1 = rising)
    
    bursts = [];
    
    dmode = diff([0 mode 0]);
    dmode(end-1) = dmode(end);
    dmode(end) = [];
    
    bursts(2,:) = find(dmode==1);
    bursts(3,:) = find(dmode==-1);
    
    if (exist('ttype', 'var'))
        bursts(5,:) = ttype(bursts(2,:));
    else
        bursts(5,:) = 0;
    end
    
    % discard the bursts that don't have any stimuli
    keeper = false(1, size(bursts, 2));
    for bursti = 1:size(bursts, 2)
        r = bursts(2,bursti):bursts(3,bursti);
        keeper(bursti) = sum(smon(1, r)) > 0;
    end
    
    bursts(:, ~keeper) = [];
    
    bursts(1,:) = 1:size(bursts,2);
    
    %% build a table with the following
    % 1 - stim id
    % 2 - sample number where occurred
    % 3 - mode
    % 4 - burst before
    % 5 - n samples after previous burst
    % 6 - burst after
    % 7 - n samples before previous burst
    % 8 - stim type (rising or falling edge of beta)
    
    stims = [];
    
    stims(2,:) = find(smon(1,:)==1);
    stims(1,:) = 1:size(stims,2);
    stims(3,:) = mode(1,stims(2,:));
    
    for stimi = 1:size(stims,2)
        if (stims(3,stimi)==1) % in burst
            stims(4:7, stimi) = NaN;
            stims(8, stimi) = ttype(stims(2,stimi));
        else
            prebursti = find(stims(2,stimi) - bursts(3,:) > 0, 1, 'last');
            if (isempty(prebursti))
                stims(4:5, stimi) = NaN;
            else
                stims(4, stimi) = prebursti;
                stims(5, stimi) = stims(2,stimi)-bursts(3,prebursti);
            end
            
            postbursti = find(bursts(2,:) - stims(2, stimi) > 0, 1, 'first');
            if (isempty(postbursti))
                stims(6:7) = NaN;
            else
                stims(6, stimi) = postbursti;
                stims(7, stimi) = bursts(2, postbursti) - stims(2, stimi);
            end
            
            stims(8, stimi) = NaN;
        end
    end
    
    %% go back to the bursts array and figure out how many ct's for each burst
    
    for bursti = 1:size(bursts, 2)
        bursts(4,bursti) = sum(stims(2,:) > bursts(2,bursti) & stims(2,:) < bursts(3, bursti));
    end
    
    % DJC - moved enumarating to AFTER, in order to account for deletions
    % DJC 9-2-2015 - get rid of any place where stim type was 2 (null
    % condition), and considered conditioning stimulation
    
    stims(:,(stims(3,:)==1 & stims(8,:)==2)) = [];
    
    stims(1,:) = 1:size(stims,2); % indexing (enumerating) each individual stim sequentially
    
    
    %% save the result to intermediate file for future use
    if saveIt
        save(fullfile(folderTiming, [sid '_tables.mat']), 'bursts', 'fs', 'stims');
    end
    
end

