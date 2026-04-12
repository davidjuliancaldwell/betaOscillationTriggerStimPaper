function T = compute_burst_phase_precision(sid)
% COMPUTE_BURST_PHASE_PRECISION  For each test stim, compute the phase
% precision of beta-band fits in the preceding conditioning burst.
%
% T = compute_burst_phase_precision('0b5a2e')
%
% Phase screening: R-squared > 0.7 AND frequency 12-20 Hz (beta band only).
% Circular mean used throughout (standard for circular data).
% nGoodRsq tracks fits with good R-squared at any frequency (diagnostic).

setup_environment;
Z_Constants;

modifierPhase = '_51samps_12_20_40ms_randomstart';

% subject config (from multipleSubj_GLMM_script_PP.m)
valueSet = {{'s',180,1,[54 62],[1 49 58 59],[44 45 46 52 53 55 60 61 63],53},...
    {'m',[0 180],2,[55 56],[1 2 3 31 57],[47 48 64],64},...
    {'s',180,3,[11 12],[57],[4 5 10 13],4},...
    {'s',270,4,[59 60],[1 9 10 35 43],[50 51 52 53 58],51},...
    {'m',[90,270],5,[13 14],[23 27 28 29 30 32 44 52 60],[5],5},...
    {'t',[270,90,12345,12345],6,[56 64],[57:63],[47 48 54 55],55},...
    {'m',[90,270],7,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31},...
    {'m',[90,270],8,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31}};
M = containers.Map(SIDS, valueSet, 'UniformValues', false);

info = M(sid);
type = info{1};
desiredF = info{2};
subjectNum = info{3};
stimElectrodes = info{4};
badChans = info{5};
goodEPs = info{6};
betaChan = info{7};

rsquaredThresh = 0.7;
freqMin = 12.01;
freqMax = 19.99;

%% load stim/burst tables
if strcmp(sid, '0b5a2ePlayback')
    load(fullfile(folderData, 'stim_timing_data', '0b5a2e_tables.mat'), 'bursts', 'fs', 'stims');
    delay = 577869;
    stims(2,:) = stims(2,:) + delay;
    bursts(2,:) = bursts(2,:) + delay;
    bursts(3,:) = bursts(3,:) + delay;
else
    load(fullfile(folderData, 'stim_timing_data', [sid '_tables.mat']), 'bursts', 'fs', 'stims');
end

stims(:, stims(2,:) < fs/2) = [];
badsStim = stims(3,:) == 0 & (isnan(stims(4,:)) | isnan(stims(6,:)));
stims(:, badsStim) = [];

%% load phase data
phaseFile = [sid '_phaseDelivery_allChans' modifierPhase '.mat'];
if strcmp(sid, '0b5a2ePlayback')
    % load the PlayBack phase file (capital B matches B_phaseCalc output)
    phaseFile = ['0b5a2ePlayBack_phaseDelivery_allChans' modifierPhase '.mat'];
end
origSid = sid;  % preserve before load (phase .mat contains 'sid' that would overwrite)
load(fullfile(folderData, 'phase_data', phaseFile));
sid = origSid;

%% build condition maps
% Maps burst types to their phase data variables.
% phase data is (trials x channels) in MATLAB.
if strcmp(type, 's')
    % d5cd55 has extra time filter in B_phaseCalc (stims(2,:) > 36536266)
    % phase_at_0 only contains fits for stims passing this filter
    condStimMask = stims(3,:) == 1;
    if strcmp(sid, 'd5cd55')
        condStimMask = condStimMask & stims(2,:) > 36536266;
    end
    condMaps = {struct('stim_idx', find(condStimMask), ...
        'phase', phase_at_0, 'rsq', r_square, 'freq', f, ...
        'burst_types', unique(bursts(5,:)), 'label', desiredF)};

elseif strcmp(type, 'm')
    % pos = stims(8)==1, neg = stims(8)==0 for c91479, 702d24, 0b5a2e
    condMaps = {
        struct('stim_idx', find(stims(3,:)==1 & stims(8,:)==1), ...
            'phase', phase_at_0_pos, 'rsq', r_square_pos, 'freq', f_pos, ...
            'burst_types', 1, 'label', desiredF(1));
        struct('stim_idx', find(stims(3,:)==1 & stims(8,:)==0), ...
            'phase', phase_at_0_neg, 'rsq', r_square_neg, 'freq', f_neg, ...
            'burst_types', 0, 'label', desiredF(2))
    };

elseif strcmp(type, 't')
    % ecb43e: 270=stims(8)==0, 90=stims(8)==1, random=stims(8)==3
    condMaps = {
        struct('stim_idx', find(stims(3,:)==1 & stims(8,:)==0), ...
            'phase', phase_at_0_pos, 'rsq', r_square_pos, 'freq', f_pos, ...
            'burst_types', 0, 'label', 270);
        struct('stim_idx', find(stims(3,:)==1 & stims(8,:)==1), ...
            'phase', phase_at_0_neg, 'rsq', r_square_neg, 'freq', f_neg, ...
            'burst_types', 1, 'label', 90);
        struct('stim_idx', find(stims(3,:)==1 & stims(8,:)==3), ...
            'phase', phase_at_0, 'rsq', r_square, 'freq', f, ...
            'burst_types', [2 3], 'label', 12345)
    };
end

%% classify probes (same logic as B_ExtractNeuralData_PP_reref.m)
% baseline: > 2s since burst end
% conditioned: < 0.5s since burst end
% excluded: 0.5-2s gap (not analyzed)
probe_mask = stims(3,:) == 0;
pstims = stims(:, probe_mask);
nProbes = size(pstims, 2);

is_baseline = pstims(5,:) > 2 * fs;
is_conditioned = pstims(5,:) < 0.5 * fs;

% dose binning: same groups as B_ExtractNeuralData_PP_reref.m
labelGroupStarts = [1 3 5];
labelGroupEnds = [labelGroupStarts(2:end) Inf];

fprintf('%s: %d probes (%d baseline, %d conditioned, %d excluded), %d channels\n', ...
    sid, nProbes, sum(is_baseline), sum(is_conditioned), ...
    sum(~is_baseline & ~is_conditioned), length(goodEPs));

%% main loop: for each (probe, channel), compute burst phase precision
results = [];

for chanIdx = 1:length(goodEPs)
    chan = goodEPs(chanIdx);

    for probeNum = 1:nProbes
        probeSample = pstims(2, probeNum);
        burstId = pstims(4, probeNum);

        if isnan(burstId) || burstId < 1 || burstId > size(bursts, 2)
            continue;
        end

        if is_baseline(probeNum)
            probeClass = 0;
        elseif is_conditioned(probeNum)
            probeClass = 1;
        else
            continue;
        end

        burstStart = bursts(2, burstId);
        burstStop = bursts(3, burstId);
        burstType = bursts(5, burstId);
        nCondInBurst = bursts(4, burstId);

        % dose label (0=Base, 1=[1,2], 2=[3,4], 3=[5,inf))
        if probeClass == 0
            doseLabel = 0;
        else
            doseLabel = 0;
            for gi = 1:length(labelGroupStarts)
                if nCondInBurst >= labelGroupStarts(gi) && nCondInBurst < labelGroupEnds(gi)
                    doseLabel = gi;
                    break;
                end
            end
        end

        % find which condition map matches this burst type
        matchedCond = 0;
        for ki = 1:length(condMaps)
            if any(burstType == condMaps{ki}.burst_types)
                matchedCond = ki;
                break;
            end
        end

        burstCircMean = NaN;
        burstVecLength = NaN;
        burstCircStd = NaN;
        nGoodRsq = 0;
        nGoodBeta = 0;

        if matchedCond > 0 && nCondInBurst > 0 && probeClass == 1
            cm = condMaps{matchedCond};

            % find conditioning stims in this burst by sample time
            condSamples = stims(2, cm.stim_idx);
            inBurst = condSamples >= burstStart & condSamples <= burstStop;
            burstCondIdx = find(inBurst);

            if ~isempty(burstCondIdx) && chan <= size(cm.phase, 2)
                % phase data is (trials x channels)
                burstPhases = cm.phase(burstCondIdx, chan)';
                burstRsq = cm.rsq(burstCondIdx, chan)';
                burstFreq = cm.freq(burstCondIdx, chan)';

                % diagnostic: how many had good R-squared at any freq
                nGoodRsq = sum(burstRsq > rsquaredThresh);

                % beta-band fits: R-squared > 0.7 AND 12-20 Hz
                betaFitMask = burstRsq > rsquaredThresh & burstFreq > freqMin & burstFreq < freqMax;
                betaPhases = burstPhases(betaFitMask);
                nGoodBeta = sum(betaFitMask);

                % circular stats on beta-band fits only
                if nGoodBeta > 0
                    burstCircMeanRad = circ_mean(betaPhases');
                    if burstCircMeanRad < 0
                        burstCircMeanRad = burstCircMeanRad + 2*pi;
                    end
                    burstCircMean = rad2deg(burstCircMeanRad);
                    burstVecLength = circ_r(betaPhases');
                    if nGoodBeta > 1
                        burstCircStd = rad2deg(circ_std(betaPhases'));
                    else
                        burstCircStd = 0;
                    end
                end
            end
        end

        targetPhase = 0;
        if matchedCond > 0
            targetPhase = condMaps{matchedCond}.label;
            if ~isscalar(targetPhase)
                targetPhase = targetPhase(1);
            end
        end

        results = [results; probeNum, probeSample, burstId, nCondInBurst, ...
            nGoodRsq, nGoodBeta, ...
            targetPhase, doseLabel, probeClass, ...
            burstCircMean, burstVecLength, burstCircStd, chan, subjectNum];
    end
end

%% build table
T = array2table(results, 'VariableNames', ...
    {'probeNum', 'probeSample', 'burstId', 'nCondStims', ...
     'nGoodRsq', 'nGoodBeta', ...
     'targetPhase', 'doseLabel', 'probeClass', ...
     'burstCircMean', 'burstVecLength', 'burstCircStd', 'channel', 'subjectNum'});

T.channelEncoded = T.subjectNum * 100 + T.channel;

doseNames = {'Base', '[1,2]', '[3,4]', '[5,inf)'};
T.doseStr = doseNames(T.doseLabel + 1)';

outfile = fullfile('data', 'output_table', [sid '_burst_phase_precision.csv']);
writetable(T, outfile);
fprintf('Saved %d rows to %s\n', height(T), outfile);

%% summary
fprintf('\n=== Summary (beta channel %d) ===\n', betaChan);
betaRows = T(T.channel == betaChan & T.probeClass == 1, :);

for targetVal = unique(betaRows.targetPhase)'
    condRows = betaRows(betaRows.targetPhase == targetVal, :);
    highDose = condRows(condRows.nCondStims >= 5, :);
    hasBeta = highDose(highDose.nGoodBeta > 0, :);

    fprintf('\ntargetPhase=%d: %d conditioned, %d at [5,inf), %d with beta fits\n', ...
        targetVal, height(condRows), height(highDose), height(hasBeta));
    fprintf('  >=1 good R-squared (any freq): %d/%d\n', sum(highDose.nGoodRsq > 0), height(highDose));
    fprintf('  >=1 beta fit (R²>0.7 & 12-20Hz): %d/%d\n', height(hasBeta), height(highDose));

    if height(hasBeta) > 0
        grandCircMeanRad = circ_mean(deg2rad(hasBeta.burstCircMean));
        grandCircMean = rad2deg(grandCircMeanRad);
        if grandCircMean < 0; grandCircMean = grandCircMean + 360; end
        grandVecLength = circ_r(deg2rad(hasBeta.burstCircMean));

        fprintf('  Beta fits per burst: mean=%.1f, median=%.0f, range=[%d,%d]\n', ...
            mean(hasBeta.nGoodBeta), median(hasBeta.nGoodBeta), ...
            min(hasBeta.nGoodBeta), max(hasBeta.nGoodBeta));
        fprintf('  Per-burst vector length: mean=%.3f, median=%.3f\n', ...
            mean(hasBeta.burstVecLength), median(hasBeta.burstVecLength));
        fprintf('  Grand circular mean: %.1f deg, grand vector length=%.3f\n', ...
            grandCircMean, grandVecLength);
    end
end

end
