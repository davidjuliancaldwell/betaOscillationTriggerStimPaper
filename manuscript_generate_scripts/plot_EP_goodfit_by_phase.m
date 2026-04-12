function plot_EP_goodfit_by_phase(sid, chanInt, varargin)
% PLOT_EP_GOODFIT_BY_PHASE  Plot average EP waveforms separated by phase
% condition and dose level, restricted to probes after good beta fits.
%
% plot_EP_goodfit_by_phase('0b5a2e', 14)
% plot_EP_goodfit_by_phase('702d24', 5, 'useMedian', true, 'saveIt', true)
%
% For each phase condition, plots dose-colored average (or median) EP
% waveforms from probes where the preceding burst had at least one
% conditioning stim with R^2 > 0.7 and frequency 12-20 Hz.
%
% Options (name-value pairs):
%   'useMedian'  - true/false (default true): use median instead of mean
%   'saveIt'     - true/false (default true): save PNG and EPS
%   'showAll'    - true/false (default true): also plot unfiltered for comparison

p = inputParser;
addRequired(p, 'sid', @ischar);
addRequired(p, 'chanInt', @isnumeric);
addParameter(p, 'useMedian', true, @islogical);
addParameter(p, 'saveIt', true, @islogical);
addParameter(p, 'showAll', true, @islogical);
parse(p, sid, chanInt, varargin{:});
opts = p.Results;

setup_environment;

%% subject config (from multipleSubj_GLMM_script_PP.m valueSet)
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
subjectNum = info{3};
desiredF = info{2};
stimElectrodes = info{4};
badChans = info{5};
betaChan = info{7};

% EP measurement window per subject (must match B_ExtractNeuralData_PP_reref.m)
tminMap = containers.Map(...
    {'d5cd55','c91479','7dbdec','9ab7ab','702d24','ecb43e','0b5a2e','0b5a2ePlayback'}, ...
    {0.006, 0.006, 0.006, 0.006, 0.00323, 0.006, 0.005, 0.005});
tmaxMap = containers.Map(...
    {'d5cd55','c91479','7dbdec','9ab7ab','702d24','ecb43e','0b5a2e','0b5a2ePlayback'}, ...
    {0.06, 0.06, 0.06, 0.06, 0.025, 0.06, 0.06, 0.06});
info_tmin = tminMap(sid);
info_tmax = tmaxMap(sid);

% peak-to-peak extraction options per subject
if strcmp(sid, '702d24')
    ppOpt = 'falling';
    ppMinPeakDist = round(0.005 * 24414);  % 5 ms
    smoothOrder = 3;
smoothFramelen = 91;

else
    ppOpt = 'abs';
    ppMinPeakDist = 15;
    smoothOrder = 3;
smoothFramelen = 171;

end

%% load stim/burst tables
if strcmp(sid, '0b5a2ePlayback')
    load(fullfile(folderData, 'stim_timing_data', '0b5a2e_tables.mat'), 'bursts', 'fs', 'stims');
    stims(2,:) = stims(2,:) + 577869;
    bursts(2,:) = bursts(2,:) + 577869;
    bursts(3,:) = bursts(3,:) + 577869;
else
    load(fullfile(folderData, 'stim_timing_data', [sid '_tables.mat']), 'bursts', 'fs', 'stims');
end
stims(:, stims(2,:) < fs/2) = [];
badsStim = stims(3,:) == 0 & (isnan(stims(4,:)) | isnan(stims(6,:)));
stims(:, badsStim) = [];
% match main pipeline's delayDelivery shift (hardware stim command → actual stim)
delayDelivery = 14;
stims(2,:) = stims(2,:) + delayDelivery;
bursts(2,:) = bursts(2,:) + delayDelivery;
bursts(3,:) = bursts(3,:) + delayDelivery;

%% load precision CSV for good-fit filtering
precFile = fullfile('data', 'output_table', [sid '_burst_phase_precision.csv']);
if ~exist(precFile, 'file')
    error('Precision CSV not found: %s. Run compute_burst_phase_precision first.', precFile);
end
precData = readtable(precFile);
% filter to this channel
chanEncoded = subjectNum * 100 + chanInt;
precChan = precData(precData.channelEncoded == chanEncoded, :);
fprintf('Loaded %d precision rows for channel %d (encoded %d)\n', ...
    height(precChan), chanInt, chanEncoded);

%% load ECoG data and re-reference
% determine re-referencing channels (all non-bad, non-stim channels)
allChans = 1:64;
badsTotal = [stimElectrodes badChans];
rerefChans = allChans(~ismember(allChans, badsTotal));

grp = floor((chanInt-1)/16);
ev = sprintf('ECO%d', grp+1);
achan = chanInt - grp*16;
load(fullfile(folderECoGData, [sid '_ECoG.mat']), ev);
dataStruct = eval(ev);
eco = 4 * dataStruct.data(:, achan)';
efs = dataStruct.info.SamplingRateHz;
fac = fs / efs;

% build re-referencing signal from median of reref channels
fprintf('Building median re-reference from %d channels...\n', length(rerefChans));
rerefData = [];
for rc = rerefChans
    grpR = floor((rc-1)/16);
    evR = sprintf('ECO%d', grpR+1);
    achanR = rc - grpR*16;
    if achanR == 1 || achanR == 2
        load(fullfile(folderECoGData, [sid '_ECoG.mat']), evR);
        dataStruct = eval(evR);
    end
    rerefData(end+1,:) = 4 * dataStruct.data(:, achanR)';
end
rerefSignal = median(rerefData, 1);
eco = eco - rerefSignal;

%% epoch around probes
probeMask = stims(3,:) == 0;
pstims = stims(:, probeMask);
probeSamples = pstims(2,:);

presamps = round(0.050 * efs);
postsamps = round(0.120 * efs);
ptis = round(probeSamples / fac);
t = (-presamps:postsamps) / efs;

wins = squeeze(getEpochSignal(eco', ptis-presamps, ptis+postsamps+1));
% baseline correct: subtract mean of pre-stim window (-50 to -5 ms)
baseWin = t < -0.005 & t > -0.050;
wins = wins - repmat(mean(wins(baseWin,:), 1), [size(wins,1), 1]);

%% classify probes: baseline vs conditioned, dose, burst type
isBaseline = pstims(5,:) > 2 * fs;
isConditioned = pstims(5,:) < 0.5 * fs;
nCondInBurst = bursts(4, pstims(4,:));
burstType = bursts(5, pstims(4,:));

% dose labels
doseLabel = zeros(1, size(pstims,2));
labelGroupStarts = [1 3 5];
labelGroupEnds = [labelGroupStarts(2:end) Inf];
for gi = 1:length(labelGroupStarts)
    mask = nCondInBurst >= labelGroupStarts(gi) & nCondInBurst < labelGroupEnds(gi) & isConditioned;
    doseLabel(mask) = gi;
end
% baseline = 0, conditioned = 1/2/3 (dose bins)

%% match probes to precision data via probeSample
% precision CSV stores probeSample BEFORE delayDelivery shift, so subtract
% delayDelivery when looking up in the precision table
goodFitMask = false(1, size(pstims,2));
for pi = 1:size(pstims,2)
    ps = probeSamples(pi) - delayDelivery;
    row = precChan(precChan.probeSample == ps, :);
    if ~isempty(row) && row.nGoodBeta(1) > 0
        goodFitMask(pi) = true;
    end
end

fprintf('Good-fit probes: %d of %d conditioned (%.0f%%)\n', ...
    sum(goodFitMask & isConditioned), sum(isConditioned), ...
    100 * sum(goodFitMask & isConditioned) / sum(isConditioned));

%% determine phase conditions
% Mirror multipleSubj_GLMM_script_PP.m: null/random burst types are skipped.
% nullType is 1-indexed into sorted unique types → bt = nullType - 1
if strcmp(sid,'0b5a2e') || strcmpi(sid,'0b5a2ePlayback') || strcmp(sid,'ecb43e')
    nullBt = 2;  % 0-indexed null burst type
else
    nullBt = NaN;
end

allTypes = unique(burstType(isConditioned));
allTypes = allTypes(~isnan(allTypes));
types = allTypes(allTypes ~= nullBt);

% Map each burst type to its desiredF target (matches
% multipleSubj_GLMM_script_PP.m:223-234 logic)
phaseLabels = cell(1, length(types));
for ti = 1:length(types)
    bt = types(ti);
    if strcmp(type,'s')
        correctIdx = 1;
    elseif strcmp(type,'m')
        correctIdx = 2 - bt;   % bt=0→2(neg), bt=1→1(pos)
    elseif strcmp(type,'t')
        if bt == 0,     correctIdx = 1;
        elseif bt == 1, correctIdx = 2;
        else,           correctIdx = 4;
        end
    end
    if correctIdx <= length(desiredF)
        phaseLabels{ti} = sprintf('%d deg', desiredF(correctIdx));
    else
        phaseLabels{ti} = sprintf('Type %d', bt);
    end
end

%% select summary function
if opts.useMedian
    avgFn = @(x, dim) median(x, dim, 'omitnan');
    statLabel = 'median';
else
    avgFn = @(x, dim) mean(x, dim, 'omitnan');
    statLabel = 'mean';
end

%% plot: individual trials (gray) + peak/trough markers + average overlay
% Layout: rows = phase x filter, columns = dose levels (4: Base, [1,2], [3,4], [5,inf))
doseNames = {'Base', '[1,2]', '[3,4]', '[5,inf)'};
doseColors = [0.5 0.5 0.5; 0.2 0.6 1.0; 0.2 0.8 0.2; 0.9 0.1 0.1];

% EP measurement window for peak-to-peak extraction
tMask = t >= info_tmin & t <= info_tmax;
tWin = t(tMask);

nPhases = length(types);
nRows = nPhases;
if opts.showAll
    nRows = nPhases * 2;
end
nCols = 4;  % Base, [1,2], [3,4], [5,inf)

fig = figure('Units', 'inches', 'Position', [0 0 4*nCols 3*nRows]);

plotIdx = 0;
for phIdx = 1:nPhases
    phaseMask = burstType == types(phIdx);

    for filterMode = 1:(1 + opts.showAll)
        if filterMode == 1
            filterLabel = 'good-fit';
            probeFilt = goodFitMask;
        else
            filterLabel = 'all';
            probeFilt = true(1, size(pstims,2));
        end

        for di = 0:3
            plotIdx = plotIdx + 1;
            subplot(nRows, nCols, plotIdx);
            hold on;

            if di == 0
                sel = isBaseline;
            else
                sel = (doseLabel == di) & phaseMask & probeFilt;
            end
            nTrials = sum(sel);
            if nTrials < 1; title(sprintf('%s (n=0)', doseNames{di+1})); continue; end

            selIdx = find(sel);
            ppVals = [];

            % plot individual trials in gray with peak/trough markers
            for ti = 1:nTrials
                trialWave = 1e6 * wins(:, selIdx(ti));
                plot(1e3*t, trialWave, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.3, ...
                    'HandleVisibility', 'off');

                % extract peak-to-peak on windowed + smoothed trial
                trialWindowed = wins(tMask, selIdx(ti));
                try
                    trialSmooth = sgolayfilt_complete(trialWindowed, smoothOrder, smoothFramelen);
                    [amp, pk_loc, tr_loc] = peak_to_peak_beta_stim(trialSmooth, ppOpt, ppMinPeakDist);
                catch
                    amp = []; pk_loc = []; tr_loc = [];
                end
                if ~isempty(amp) && ~isnan(amp)
                    % peak marker: red downward triangle
                    plot(1e3*tWin(pk_loc), 1e6*trialSmooth(pk_loc), 'v', ...
                        'Color', [0.9 0.1 0.1], 'MarkerSize', 4, ...
                        'MarkerFaceColor', [0.9 0.1 0.1], 'HandleVisibility', 'off');
                    % trough marker: blue upward triangle
                    plot(1e3*tWin(tr_loc), 1e6*trialSmooth(tr_loc), '^', ...
                        'Color', [0.1 0.1 0.9], 'MarkerSize', 4, ...
                        'MarkerFaceColor', [0.1 0.1 0.9], 'HandleVisibility', 'off');
                    ppVals(end+1) = amp * 1e6;
                end
            end

            % overlay average in color
            waveAvg = avgFn(1e6 * wins(:, sel), 2);
            plot(1e3*t, waveAvg, 'Color', doseColors(di+1,:), 'LineWidth', 2, ...
                'DisplayName', sprintf('%s (n=%d)', doseNames{di+1}, nTrials));

            xlim([-5 60]);
            ylim([-400 400]);
            if di == 0
                ylabel(sprintf('%s\n[%s]\n\\muV', phaseLabels{phIdx}, filterLabel));
            end
            if ~isempty(ppVals)
                title(sprintf('%s (n=%d)\nPP: %.0f \\muV', doseNames{di+1}, nTrials, ...
                    median(ppVals)), 'FontSize', 9);
            else
                title(sprintf('%s (n=%d)', doseNames{di+1}, nTrials), 'FontSize', 9);
            end
            set(gca, 'FontSize', 8);
        end
    end
end
sgtitle(sprintf('%s ch%d — individual trials + %s EP with peak/trough markers', ...
    sid, chanInt, statLabel), 'FontSize', 13);

%% save
if opts.saveIt
    outBase = fullfile(folderPlots, sprintf('EP_goodfit_%s_ch%d_%s', sid, chanInt, statLabel));
    saveas(fig, [outBase '.png']);
    print(fig, [outBase '.eps'], '-depsc', '-r600');
    fprintf('Saved: %s.{png,eps}\n', outBase);
end

end
