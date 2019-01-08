%% Constants
%close all; clear all;clc

Z_Constants;
SUB_DIR = fullfile(myGetenv('subject_dir'));

%% Load in the trigger data
%DJC 7/20/2015 - changed tp to fit David paths

for index = 1:8
    % select the subject from list
    sid = SIDS{index};
    
    if (strcmp(sid, 'd5cd55'))
        
        tp = strcat(SUB_DIR,'\d5cd55\data\D8\d5cd55_BetaTriggeredStim');
        block = 'Block-49';
        
        %load([block '.mat'],Smon
        % SMon-1 is the system enable
        % SMon-2 is the stim command
        % SMon-3 is the stim count
        % SMon-4 is the realized voltage
        tic;
        [smon, fs] = tdt_loadStream(tp, block, 'SMon', 2);
        [stim, ~] = tdt_loadStream(tp, block, 'SMon', 4);
        toc;
        
        % in the recording for d5cd55
        % Wave-1 looks like the trigger signal
        % Wave-2 is the mode
        % Wave-3 is the mode time/counter
        % Wave-4 looks like the stim command
        tic;
        [mode, ~] = tdt_loadStream(tp, block, 'Wave', 2);
        [beta, ~] = tdt_loadStream(tp, block, 'Wave', 1);
        ttype = 0*beta;
        toc;
        
        %     % these three lines will get rid of the large simuli at the beginning
        %     % of the record
        %     mode(1:4.5e6) = [];
        %     beta(1:4.5e6) = [];
        %     smon(1:4.5e6) = [];
        
        %     % these three lines will get rid of all stimuli until we reset the
        %     % threshold value
        %     mode(1:36536266) = [];
        %     beta(1:36536266) = [];
        %     smon(1:36536266) = [];
        
        %     % these lines will get rid of all stimuli after we reset the threshold
        %     % value
        %     mode(36536266:end) = [];
        %     beta(36536266:end) = [];
        %     smon(36536266:end) = [];
        %     mode(1:4.5e6) = [];
        %     beta(1:4.5e6) = [];
        %     smon(1:4.5e6) = [];
        
    elseif (strcmp(sid, 'c91479'))
        tp = strcat(SUB_DIR,'\c91479\data\d7\c91479_BetaTriggeredStim');
        
        block = 'BetaPhase-14';
        
        % SMon-1 is the system enable
        % SMon-2 is the stim command
        % SMon-3 is the stim count
        % SMon-4 is the realized voltage
        tic;
        [smon, fs] = tdt_loadStream(tp, block, 'SMon', 2);
        [stim, ~] = tdt_loadStream(tp, block, 'SMon', 4);
        toc;
        
        % in the recording for c91479
        % Wave-1 looks like phase decision variable (0=falling, 1=rising)
        % Wave-2 is the mode
        % Wave-3 is the mode time/counter
        % Wave-4 looks like the stim command
        tic;
        [mode, ~] = tdt_loadStream(tp, block, 'Wave', 2);
        [ttype, ~] = tdt_loadStream(tp, block, 'Wave', 1);
        [beta, ~] = tdt_loadStream(tp, block, 'Blck', 1);
        [raw, ~] = tdt_loadStream(tp, block, 'Blck', 2);
        toc;
        
        % these lines will get rid of the time period at the end of the record
        % where we were trying a lower max stim frequency
        mode(64507402:end) = 0;
        ttype(64507402:end) = 0;
        beta(64507402:end) = 0;
        smon(64507402:end) = 0;
        stim(64507402:end) = 0;
        raw(64507402:end) = 0;
        
        % these lines will get rid of the time period at the beginning of
        % the record where we were changing parameter settings
        mode(1:2e7) = 0;
        ttype(1:2e7) = 0;
        beta(1:2e7) = 0;
        smon(1:2e7) = 0;
        stim(1:2e7) = 0;
        raw(1:2e7) = 0;
    elseif (strcmp(sid, '7dbdec'))
        tp = strcat(SUB_DIR,'\7dbdec\data\d7\7dbdec_BetaTriggeredStim');
        
        block = 'BetaPhase-17';
        
        % SMon-1 is the system enable
        % SMon-2 is the stim command
        % SMon-3 is the stim count
        % SMon-4 is the realized voltage
        tic;
        [smon, fs] = tdt_loadStream(tp, block, 'SMon', 2);
        [stim, ~] = tdt_loadStream(tp, block, 'SMon', 4);
        toc;
        
        % Wave-1 looks like phase decision variable (0=falling, 1=rising)
        % Wave-2 is the mode
        % Wave-3 is the mode time/counter
        % Wave-4 looks like the stim command
        tic;
        [mode, ~] = tdt_loadStream(tp, block, 'Wave', 2);
        [ttype, ~] = tdt_loadStream(tp, block, 'Wave', 1);
        [beta, ~] = tdt_loadStream(tp, block, 'Blck', 1);
        [raw, ~] = tdt_loadStream(tp, block, 'Blck', 2);
        toc;
    elseif (strcmp(sid, '9ab7ab'))
        tp = strcat(SUB_DIR,'\9ab7ab\data\d7\9ab7ab_BetaTriggeredStim');
        
        block = 'BetaPhase-3';
        
        % SMon-1 is the system enable
        % SMon-2 is the stim command
        % SMon-3 is the stim count
        % SMon-4 is the realized voltage
        tic;
        [smon, fs] = tdt_loadStream(tp, block, 'SMon', 2);
        [stim, ~] = tdt_loadStream(tp, block, 'SMon', 4);
        toc;
        
        % Wave-1 looks like the beta signal
        
        % Wave-2 is the mode
        % Wave-3 is the mode time/counter
        % Wave-4 looks like the stim command
        tic;
        [mode, ~] = tdt_loadStream(tp, block, 'Wave', 2);
        %     [ttype, ~] = tdt_loadStream(tp, block, 'Wave', 1);
        
        [x1, ~] = tdt_loadStream(tp, block, 'Wave', 3);
        [x2, ~] = tdt_loadStream(tp, block, 'Wave', 4);
        
        ttype = 0*mode;
        [beta, ~] = tdt_loadStream(tp, block, 'Blck', 1);
        [raw, ~] = tdt_loadStream(tp, block, 'Blck', 2);
        toc;
    elseif (strcmp(sid, '702d24'))
        tp = strcat(SUB_DIR,'\702d24\data\d7\702d24_BetaStim');
        
        tank = TTank;
        tank.openTank(tp);
        tank.selectBlock('BetaPhase-4');
        %     tp = 'd:\research\subjects\702d24\data\d7\702d24_BetaStim';
        %     block = 'BetaPhase-4';
        
        % SMon-1 is the system enable
        % SMon-2 is the stim command
        % SMon-3 is the stim count
        % SMon-4 is the realized voltage
        tic;
        [smon, info] = tank.readWaveEvent('SMon', 2);
        smon = smon';
        
        fs = info.SamplingRateHz;
        
        stim = tank.readWaveEvent('SMon', 4)';
        toc;
        
        tic;
        mode = tank.readWaveEvent('Wave', 2)';
        ttype = tank.readWaveEvent('Wave', 1)';
        
        beta = tank.readWaveEvent('Blck', 1)';
        %     [beta, ~] = tdt_loadStream(tp, block, 'Blck', 1);
        raw = tank.readWaveEvent('Blck', 2)';
        %     [raw, ~] = tdt_loadStream(tp, block, 'Blck', 2);
        toc;
        
        % added last subject, DJC - 7-23-2015
    elseif (strcmp(sid, 'ecb43e'))
        tp = strcat(SUB_DIR,'\ecb43e\data\d7\BetaStim');
        
        tank = TTank;
        tank.openTank(tp);
        tank.selectBlock('BetaPhase-3');
        
        tic;
        [smon, info] = tank.readWaveEvent('SMon', 2);
        smon = smon';
        
        fs = info.SamplingRateHz;
        
        stim = tank.readWaveEvent('SMon', 4)';
        toc;
        
        % Wave-2 is the mode
        % Wave-3 is the mode time/counter
        % Wave-4 looks like the stim command
        tic;
        mode = tank.readWaveEvent('Wave', 2)';
        ttype = tank.readWaveEvent('Wave', 1)';
        
        beta = tank.readWaveEvent('Blck', 1)';
        %     [beta, ~] = tdt_loadStream(tp, block, 'Blck', 1);
        raw = tank.readWaveEvent('Blck', 2)';
        %     [raw, ~] = tdt_loadStream(tp, block, 'Blck', 2);
        toc;
        
        %     toc;
    elseif (strcmp(sid, '0b5a2e'))
        tp = strcat(SUB_DIR,'\0b5a2e\data\d8\0b5a2e_BetaStim\0b5a2e_BetaStim');
        
        tank = TTank;
        tank.openTank(tp);
        tank.selectBlock('BetaPhase-2');
        
        tic;
        [smon, info] = tank.readWaveEvent('SMon', 2);
        smon = smon';
        
        fs = info.SamplingRateHz;
        
        stim = tank.readWaveEvent('SMon', 4)';
        toc;
        
        tic;
        mode = tank.readWaveEvent('Wave', 2)';
        ttype = tank.readWaveEvent('Wave', 1)';
        
        beta = tank.readWaveEvent('Blck', 1)';
        
        raw = tank.readWaveEvent('Blck', 2)';
        
        toc;
    elseif (strcmp(sid, '0b5a2ePlayback'))
        tp = strcat(SUB_DIR,'\0b5a2e\data\d8\0b5a2e_BetaStim\0b5a2e_BetaStim');
        
        tank = TTank;
        tank.openTank(tp);
        
        tank.selectBlock('BetaPhase-4');
        
        tic;
        [smon, info] = tank.readWaveEvent('SMon', 2);
        smon = smon';
        
        fs = info.SamplingRateHz;
        
        stim = tank.readWaveEvent('SMon', 4)';
        toc;
        
        tic;
        mode = tank.readWaveEvent('Wave', 2)';
        ttype = tank.readWaveEvent('Wave', 1)';
        
        beta = tank.readWaveEvent('Blck', 1)';
        
        raw = tank.readWaveEvent('Blck', 2)';
        
        toc;
        
    else
        error('unknown sid entered');
    end
    
    save([sid '_stim_table_data_file'],'smon','fs','stim','mode','ttype','beta')
end