%% script to plot the peak to peak amplitude differences vs. actual phase of delivery
%
% David.J.Caldwell 10.2.2018

setup_environment

SIDS = {'d5cd55','c91479','7dbdec','9ab7ab','702d24','ecb43e','0b5a2e','0b5a2ePlayback'};
% valueSet = {{'s',180,1,[54 62],[1 49 58 59],[44 45 46 47 48 52 53 55 60 61 63],53},...
%     {'m',[0 180],2,[55 56],[1 2 3 31 57],[39 40 47 48 63 64],64},...
%     {'s',180,3,[11 12],[57],[4 5 10 13 18 19 20],4},...
%     {'s',270,4,[59 60],[1 9 10 35 43],[41 42 43 44 45 49 50 51 52 53 57 58 61 62],51},...
%     {'m',[90,270],5,[13 14],[23 27 28 29 30 32 44 52 60],[5],5},...
%     {'t',[90,270],6,[56 64],[57:64],[46 48 54 55 63],55},...
%     {'m',[90,270],7,[22 30],[24 25 29],[13 14 15 16 20 21 23 31 32 39 40],31},...
%     {'m',[90,270],8,[22 30],[24 25 29],[13 14 15 16 20 21 23 31 32 39 40],31}};

valueSet = {{'s',180,1,[54 62],[1 49 58 59],[44 45 46 52 53 55 60 61 63],53,2.5},...
    {'m',[0 180],2,[55 56],[1 2 3 31 57],[47 48 64],64,3},...
    {'s',180,3,[11 12],[57],[4 5 10 13],4,3.5},...
    {'s',270,4,[59 60],[1 9 10 35 43],[50 51 52 53 58],51,0.75},...
    {'m',[90,270],5,[13 14],[23 27 28 29 30 32 44 52 60],[5],5,0.75},...
    {'t',[270,90,12345,12345],6,[56 64],[57:63],[47 48 54 55],55,1.75}...
    {'m',[90,270],7,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31,1.75},...
    {'m',[90,270],8,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31,1.75}};

M = containers.Map(SIDS,valueSet,'UniformValues',false);
plotColor = [
    [.65, .65, .65];...   % light gray         (0)
    [0.1, 0.74, 0.95];...  % deep sky-blue     (1)
    [0.95, 0.88, 0.05];... % gold/yellow       (2)
    [0.80, 0.05, 0.78];... % magenta           (3)
    [0.3, 0.8, 0.20];...   % lime green        (4)
    [0.95, 0.1, 0.1];...   % crimson red       (5)
    [0.64, 0.18, 0.93];... % blue-violet       (6)
    [0.88, 0.56, 0];...    % orange            (7)
    [0.4, 1.0, 0.7];...    % aquamarine        (8)
    [0.95, 0.88, 0.7];...  % salmon-yellow     (9)
    [0, 0.2, 1];...        % blue              (10)
    [1, 0.41, 0.7];...     % hot pink          (11)
    [0.5, 1, 0];...        % chartreuse        (12)
    [0.6, 0.39, 0.8];...   % amtheyist         (13)
    [0.82, 0.36, 0.36,];...% indian red        (14)
    [0.53, 0.8, 0.98];...  % light sky blue    (15)
    [0, 0.6, 0.1];...      % forest green      (16)
    [0.65, 0.95, 0.5];...  % light green       (17)
    [0.85, 0.6, 0.88];...  % light purple      (18)
    [0.90, 0.7, 0.7];...   % light red         (19)
    [0.2, 0.2, 0.6];...    % dark blue         (20)
    ];

% djc 4/8/26 - including 702d24 again
SIDSint = {'d5cd55','c91479','7dbdec','9ab7ab','702d24','ecb43e','0b5a2e'};

plotColor = distinguishable_colors(9);

%modifierPhase = '_13samps_8_30_40ms_randomstart';

modifierPhase = '_51samps_12_20_40ms_randomstart';

%modifierPhase = '_13samps_10_30_40ms_randomstart';

modifierEP = '-reref-50-new';
%SIDS = {'d5cd55'};

% decide how to plot circles - std deviation or vector length
markerToUse = 'vecLength';
testStatistic = 'omnibus';

threshold = 0.7;
fThresholdMin = 12.01;
fThresholdMax = 19.99;
%
% fThresholdMin = 10;
% fThresholdMax = 29.99;
markerMin = 50;
markerMax = 500;
minData = 0;
maxData = 1;
epThresholdMaxMean = 100;
epThresholdMin = 25;
epThresholdMax = 1500;

% --- Phase quality / good-fit filters (configurable; mirror R Model 5a / 5a-gf) ---
% If these variables already exist in the workspace (e.g., set by a
% wrapper script before calling `run('phase_vs_peak.m')`), the existing
% values are respected. Otherwise the defaults below apply.
%
% To generate both Model 5a and Model 5a-gf2 variants in one session, run
% the script twice via generate_both_phase_vs_peak.m or by pre-setting
% these variables in the MATLAB workspace before calling `run()`.
%
% R defaults: 0.3 (Model 5a), 0.2 (Model 5a-gf2).
if ~exist('minPhaseVecLength', 'var'),   minPhaseVecLength = 0; end
if ~exist('applyGoodFitFilter', 'var'),  applyGoodFitFilter = false; end
if ~exist('minGoodBetaPerBurst', 'var'), minGoodBetaPerBurst = 1; end
if ~exist('minBurstVecLength', 'var'),   minBurstVecLength = 0; end
if ~exist('plotExcluded', 'var'),        plotExcluded = false; end

% summary statistic: 'mean' or 'median'
% median is consistent with the summary-level LME pipeline (Models 3-5)
summaryStatType = 'median';
if strcmp(summaryStatType, 'median')
    avgFn = @(x) nanmedian(x);
else
    avgFn = @(x) nanmean(x);
end


%% plot EP modulation vs phase for all subjects
figTotal = figure;
hold on
figInd = figure;
countInd = 1;
h = gobjects(1, length(SIDSint));
hold on
for sid = SIDSint
    
    sid = sid{:};
    subjid = sid;
    info = M(sid);
    type = info{1};
    subjectNum = info{3};
    desiredF = info{2};
    stims = info{4};
    bads = info{5};
    goodEPs = info{6};
    betaChan = info{7};
    chans = [1:64];
    badsTotal = [stims bads];
    chans(ismember(chans, badsTotal) | ~ismember(chans,goodEPs)) = [];
    Montage.MontageTokenized = {'Grid(1:64)'};
    
    wInd = [];
    %  h = [];
    peakPhaseVec = [];
    chanVec = [];
    peakPhaseRep = [];
    indexVec = [];
    
    load(fullfile(folderEP,strcat(subjid,['epSTATS-PP-sig' modifierEP '.mat'])))
    load(fullfile(folderPhase,[sid '_phaseDelivery_allChans' modifierPhase '.mat']));

    fprintf(['running for subject ' sid '\n']);

    % --- Optional: load stim table and precision CSV for good-fit trial filter ---
    % Probe ordering in dataForPPanalysis{chan}{ii}{1} matches the pts selector
    % in B_ExtractNeuralData_PP_reref.m, which (1) drops first-500ms stims,
    % (2) drops bad probes, (3) applies delayDelivery=+14 sample shift,
    % (4) selects stims(3,:)==0 — plus, for d5cd55 only, an extra
    % shifted_sample > 36536266 filter. We must replicate that pts selector
    % here (in un-shifted coordinates, matching compute_burst_phase_precision.m)
    % so probeSampleTimesSID(ti) is the correct probe sample for trial ti.
    if applyGoodFitFilter
        if strcmp(sid, '0b5a2ePlayback')
            stimTable = load(fullfile(folderTiming,'0b5a2e_tables.mat'),'stims','fs');
            stimTable.stims(2,:) = stimTable.stims(2,:) + 577869;
        else
            stimTable = load(fullfile(folderTiming,[sid '_tables.mat']),'stims','fs');
        end
        stimTable.stims(:, stimTable.stims(2,:) < stimTable.fs/2) = [];
        badProbes = stimTable.stims(3,:)==0 & ...
            (isnan(stimTable.stims(4,:)) | isnan(stimTable.stims(6,:)));
        stimTable.stims(:, badProbes) = [];

        delayDelivery_extract = 14;
        if strcmp(sid, 'd5cd55')
            d5cd55ProbeThresh = 36536266 - delayDelivery_extract;  % = 36536252 un-shifted
            pts_extract = stimTable.stims(3,:) == 0 & stimTable.stims(2,:) > d5cd55ProbeThresh;
        else
            pts_extract = stimTable.stims(3,:) == 0;
        end
        probeStimsLocal = stimTable.stims(:, pts_extract);
        probeSampleTimesSID = probeStimsLocal(2,:);  % un-shifted samples

        precFile = fullfile(folderOutput,[sid '_burst_phase_precision.csv']);
        if exist(precFile,'file')
            precData = readtable(precFile);
        else
            warning('Precision CSV missing for %s (looked for %s) — good-fit filter cannot be applied; channels will be excluded from filtered plot', ...
                sid, precFile);
            precData = [];
        end
    else
        probeSampleTimesSID = [];
        precData = [];
    end

    %%

    if strcmp(type,'m')
        indices = [1,2];
    elseif strcmp(type,'s')
        indices = 1;
    elseif strcmp(type,'t')
        indices = [1,2,4];
    end
    w = nan(length(chans), length(indices));
    wExcl = nan(length(chans), length(indices));  % excluded channels (unfiltered y) for open markers
    
    for index = indices

        % Burst-type ↔ phase variable pairing.
        %
        % B_ExtractNeuralData_PP_reref.m indexes dataForPPanalysis{chan}{typei}
        % by typei = 1:length(types) where types = sort(unique(bursts(5,:))).
        %
        % For 'm' type (c91479, 702d24, 0b5a2e) the hardware convention is
        % stims(8)==0 → neg and stims(8)==1 → pos, so sorted types = [0, 1]
        % gives typei=1 → burst type 0 (neg, target desiredF(2)) and
        % typei=2 → burst type 1 (pos, target desiredF(1)).
        %
        % For 't' type (ecb43e) the convention is INVERTED: stims(8)==0 →
        % pos(270°), stims(8)==1 → neg(90°). Sorted types = [0,1,2,3] still
        % gives typei=1 → burst type 0 and typei=2 → burst type 1, but now
        % typei=1 is the POS side and typei=2 is the NEG side. typei=3 is
        % the null burst (skipped). typei=4 is the random condition.
        %
        % This is the same 2026-04-06 fix that was applied to
        % multipleSubj_GLMM_script_PP.m via `correctIdx = 2 - bt`. It was
        % not backported to phase_vs_peak.m until 2026-04-11 — before that,
        % the 'm' type branches below had phase_at_0_pos/neg SWAPPED
        % against the extraction's stored burst types, so the y-values of
        % the two dots per channel (for multi-phase subjects) were paired
        % with the wrong x-coordinates.
        if strcmp(type, 's')
            rsq_use   = r_square;
            phase_use = phase_at_0;
            f_use     = f;
            target    = desiredF;
        elseif strcmp(type, 'm') && index == 1
            % typei=1 = burst type 0 (neg), target desiredF(2)
            rsq_use   = r_square_neg;
            phase_use = phase_at_0_neg;
            f_use     = f_neg;
            target    = desiredF(2);
        elseif strcmp(type, 'm') && index == 2
            % typei=2 = burst type 1 (pos), target desiredF(1)
            rsq_use   = r_square_pos;
            phase_use = phase_at_0_pos;
            f_use     = f_pos;
            target    = desiredF(1);
        elseif strcmp(type, 't') && index == 1
            % typei=1 = burst type 0 (pos under inverted convention), target desiredF(1)=270°
            rsq_use   = r_square_pos;
            phase_use = phase_at_0_pos;
            f_use     = f_pos;
            target    = desiredF(1);
        elseif strcmp(type, 't') && index == 2
            % typei=2 = burst type 1 (neg under inverted convention), target desiredF(2)=90°
            rsq_use   = r_square_neg;
            phase_use = phase_at_0_neg;
            f_use     = f_neg;
            target    = desiredF(2);
        elseif strcmp(type, 't') && index == 4
            % typei=4 = burst type 3 = random (stims(8)==3); target = 12345
            rsq_use   = r_square;
            phase_use = phase_at_0;
            f_use     = f;
            target    = 12345;
        else
            error('phase_vs_peak: no phase mapping for type=%s index=%d', type, index);
        end

        [peakPhase, peakStd, peakLength, circularTest, markerSize] = ...
            phase_delivery_accuracy_forPP(rsq_use, threshold, phase_use, chans, ...
                target, markerMin, markerMax, minData, maxData, markerToUse, ...
                testStatistic, f_use, fThresholdMin, fThresholdMax);

        peakPhaseVec(index,:) = peakPhase;
        
        count = 1;
        for i = chans
            mags = 1e6*dataForPPanalysis{i}{index}{1};
            mags(mags<epThresholdMin) = nan;
            mags(mags>epThresholdMax) = nan;

            label = dataForPPanalysis{i}{index}{4};
            keeps = dataForPPanalysis{i}{index}{5};
            maxLabel = max(unique(label)); % plot vs. maximum number of stimuli tested
            baseVal = avgFn(mags(label ==0 & keeps));

            % unfiltered condition value (all probes passing keeps + mag bounds)
            condMaskRaw = (label == maxLabel) & keeps;
            condValRaw = avgFn(mags(condMaskRaw));
            diffRaw = 100*(condValRaw - baseVal)/baseVal;

            % optionally apply trial-level good-fit filter (Model 5a-gf analog)
            if applyGoodFitFilter
                if isempty(precData)
                    % Precision CSV missing for this subject (already warned
                    % at load time). Mark NaN so channel is excluded from the
                    % filtered plot rather than silently using unfiltered vals.
                    diffFilt = NaN;
                else
                    chanEnc = subjectNum*100 + i;
                    precChanRows = precData(precData.channelEncoded == chanEnc, :);
                    if isempty(precChanRows)
                        % Channel not processed by compute_burst_phase_precision
                        % (e.g. not in its goodEPs list). Fail loudly instead
                        % of silently pretending the filter was applied.
                        warning('phase_vs_peak:noPrecRows', ...
                            '%s ch%d: no rows in %s_burst_phase_precision.csv — regenerate precision CSV or update goodEPs. Channel excluded from filtered plot.', ...
                            sid, i, sid);
                        diffFilt = NaN;
                    else
                        % probeSample -> (nGoodBeta, burstVecLength) maps
                        keysVec = num2cell(double(precChanRows.probeSample));
                        ngMap   = containers.Map(keysVec, num2cell(double(precChanRows.nGoodBeta)));
                        bvlMap  = containers.Map(keysVec, num2cell(double(precChanRows.burstVecLength)));
                        condMaskFilt = false(size(condMaskRaw));
                        condIdx = find(condMaskRaw);
                        nMissingLookup = 0;
                        for ti = condIdx(:)'
                            if ti > length(probeSampleTimesSID)
                                % extraction has more trials than probe vec;
                                % should never happen after the d5cd55 fix
                                nMissingLookup = nMissingLookup + 1;
                                continue;
                            end
                            psLookup = double(probeSampleTimesSID(ti));
                            if ~isKey(ngMap, psLookup)
                                % probe exists in extraction but not precision
                                % CSV — alignment drift or stale CSV
                                nMissingLookup = nMissingLookup + 1;
                                continue;
                            end
                            if ngMap(psLookup) < minGoodBetaPerBurst; continue; end
                            if minBurstVecLength > 0 && bvlMap(psLookup) < minBurstVecLength
                                continue;
                            end
                            condMaskFilt(ti) = true;
                        end
                        if nMissingLookup > 0
                            warning('phase_vs_peak:probeLookupMiss', ...
                                '%s ch%d: %d of %d conditioned trials had no matching precision row (alignment drift or stale CSV).', ...
                                sid, i, nMissingLookup, length(condIdx));
                        end
                        condValFilt = avgFn(mags(condMaskFilt));
                        diffFilt = 100*(condValFilt - baseVal)/baseVal;
                    end
                end
            else
                diffFilt = diffRaw;
            end

            % channel-level phase vector length filter (Model 5a analog)
            chanLenOk = (minPhaseVecLength <= 0) || ...
                (~isnan(peakLength(count)) && peakLength(count) >= minPhaseVecLength);
            baselineOk = baseVal > epThresholdMaxMean;

            if baselineOk && chanLenOk
                w(count,index) = diffFilt;
                wTotal(subjectNum,count,index) = diffFilt;
                phaseTotal(subjectNum,count,index) = peakPhase(count);
            else
                w(count,index) = nan;
                wTotal(subjectNum,count,index) = nan;
                phaseTotal(subjectNum,count,index) = nan;
                % remember the channel's unfiltered value so it can still be
                % drawn as an open marker when plotExcluded is true
                if baselineOk
                    wExcl(count,index) = diffRaw;
                end
            end
            count = count + 1;
        end
        
        %         peakPhaseRep = repmat(peakPhaseVec',1,1,size(wInd,3));
        %
        %         chanVec = repmat([1:size(wInd,1)]',1,size(wInd,2),size(wInd,3));
        %         indexVec = repmat([1:size(wInd,2)]',size(wInd,1),1,size(wInd,3));
        %
        %         dataTable = table(wInd(:),peakPhaseRep(:),chanVec(:),indexVec(:));
        %         dataFit = fitlm(dataTable,'Var1~Var2');
        %
        
        
        %         hold on
        %         h = plot(dataFit);
        %         s=findobj('type','legend');
        %         h(1).Color = plotColor(subjectNum,:);
        %         h(1).Marker = 'o';
        %         h(1).MarkerFaceColor = plotColor(subjectNum,:);
        %         delete(s)
        %         xlabel('');
        %         ylabel('');
        %         xlim([0 360])
        %         ylim([-30 60])
        %         xticks([0 45 90 135 180 225 270 315 360])
        %         title(['Subject '  num2str(subjectNum)])
        %         set(gca,'fontsize',14)
        
        % --- scatter: filled = included; open = excluded (if plotExcluded) ---
        includedMask = ~isnan(w(:,index));
        excludedMask = ~isnan(wExcl(:,index));

        figure(figTotal)
        hold on
        if any(includedMask)
            h(countInd) = scatter(peakPhase(includedMask), w(includedMask,index), ...
                markerSize(includedMask), plotColor(subjectNum,:), 'filled');
        else
            % placeholder so legend handle stays valid for this subject
            h(countInd) = scatter(nan, nan, markerMin, plotColor(subjectNum,:), 'filled');
        end
        if plotExcluded && any(excludedMask)
            scatter(peakPhase(excludedMask), wExcl(excludedMask,index), ...
                markerSize(excludedMask), plotColor(subjectNum,:), ...
                'LineWidth', 1.2, 'HandleVisibility','off');
        end

        figure(figInd)
        grid on
        hold on
        subplot(4,2,countInd)
        hold on
        ylim([-30 60])
        if any(includedMask)
            scatter(peakPhase(includedMask), w(includedMask,index), ...
                markerSize(includedMask), plotColor(subjectNum,:), 'filled');
        end
        if plotExcluded && any(excludedMask)
            scatter(peakPhase(excludedMask), wExcl(excludedMask,index), ...
                markerSize(excludedMask), plotColor(subjectNum,:), 'LineWidth', 1.2);
        end
        xlim([0 360])
        ylim([-30 60])
        hline(0,'k')

        xticks([0 45 90 135 180 225 270 315 360])
        title(['Subject '  num2str(subjectNum)])
        set(gca,'fontsize',14)
    end
    
    countInd = countInd + 1;
    
end
%%
figure(figTotal)
grid on
xlim([0 360])
xticks([0 45 90 135 180 225 270 315 360])
ylim([-30 60])
hline(0,'k')
legend(h,{'Subject 1',...
    'Subject 2',...
    'Subject 3',...
    'Subject 4',...
    'Subject 5',...
    'Subject 6',...
    'Subject 7'})
title('Phase of delivery and CEP modulation')
xlabel('Phase of delivery (degrees)')
ylabel({'EP percent change from baseline','to >5 conditioning stimuli'})
set(gca,'fontsize',24)

% add in scale dots for markers (circle size = vector length, range 0 to 1)
scatter(270,53,markerMin,[0.4 0.4 0.4],'filled','HandleVisibility','off')
text(278,53,'VL = 0','fontsize',14,'VerticalAlignment','middle')
scatter(270,43,markerMax,[0.4 0.4 0.4],'filled','HandleVisibility','off')
text(278,43,'VL = 1','fontsize',14,'VerticalAlignment','middle')

figure(figInd)
xlabel('Phase of delivery (degrees)')
ylabel([{'EP percent change from baseline',' to >5 conditioning stimuli'}])

% save figures
% build filter suffix so different filter settings don't overwrite each other
filterSuffix = '';
if minPhaseVecLength > 0
    filterSuffix = [filterSuffix sprintf('_r%02d', round(100*minPhaseVecLength))];
end
if applyGoodFitFilter
    filterSuffix = [filterSuffix sprintf('_gf%d', minGoodBetaPerBurst)];
    if minBurstVecLength > 0
        filterSuffix = [filterSuffix sprintf('_bvl%02d', round(100*minBurstVecLength))];
    end
end

figure(figTotal)
saveas(figTotal, fullfile(folderPlots, sprintf('phase_vs_peak_all_subj_%s%s.png', summaryStatType, filterSuffix)));
print(figTotal, fullfile(folderPlots, sprintf('phase_vs_peak_all_subj_%s%s.eps', summaryStatType, filterSuffix)), '-depsc', '-r600');

figure(figInd)
saveas(figInd, fullfile(folderPlots, sprintf('phase_vs_peak_per_subj_%s%s.png', summaryStatType, filterSuffix)));
print(figInd, fullfile(folderPlots, sprintf('phase_vs_peak_per_subj_%s%s.eps', summaryStatType, filterSuffix)), '-depsc', '-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the difference between 0-180 and 180-360
wTotal = wTotal(1:7,:,:);
phaseTotal = phaseTotal(1:7,:,:);

phaseTotalLess = phaseTotal((phaseTotal < 180) & (phaseTotal>0));
phaseTotalMore = phaseTotal((phaseTotal > 180) & (phaseTotal<365) );
wTotalLess = wTotal((phaseTotal < 180) & (phaseTotal>0));
wTotalMore = wTotal((phaseTotal > 180) & (phaseTotal<365) );


wTotalLess = wTotalLess(~isnan(wTotalLess));
wTotalMore = wTotalMore(~isnan(wTotalMore));
[h,p] = ttest2(wTotalLess,wTotalMore)

[p,h,stats] = ranksum(wTotalLess,wTotalMore)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE (2026-04-11): all sections BELOW the return above are exploratory
% / dead code. They still contain the same burst-type ↔ phase variable
% swap that was fixed in the main loop at line ~208 on 2026-04-11. Any
% 'm' type subject (c91479, 702d24, 0b5a2e, 0b5a2ePlayback) plotted by
% the sections below will have the y-values of its two dots per channel
% paired with the wrong measured phases. If you ever reactivate one of
% these sections (by removing the `return` or running a block
% interactively), apply the same `correctIdx = 2 - bt` fix. See the
% commented pairing block in the main loop above (lines ~208+) for the
% correct phase_at_0_pos/neg mapping per burst type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot EP modulation vs phase for subj. 7 with playback

figure
clearvars hPlayback1
hold on
countScatter = 1;

for sid = SIDS(end-1:end)
    sid = sid{:};
    subjid = sid;
    info = M(sid);
    type = info{1};
    subjectNum = info{3};
    desiredF = info{2};
    stims = info{4};
    bads = info{5};
    goodEPs = info{6};
    betaChan = info{7};
    chans = [1:64];
    badsTotal = [stims bads];
    chans(ismember(chans, badsTotal) | ~ismember(chans,goodEPs)) = [];
    Montage.MontageTokenized = {'Grid(1:64)'};
    
        
    load(fullfile(folderEP,strcat(subjid,['epSTATS-PP-sig' modifierEP '.mat'])))
    load(fullfile(folderPhase,[sid '_phaseDelivery_allChans' modifierPhase '.mat']));
    fprintf(['running for subject ' sid '\n']);
    
    %%
    if strcmp(type,'m')
        indices = [1,2];
    elseif strcmp(type,'s')
        indices = 1;
    elseif strcmp(type,'t')
        indices = [1,2,4];
    end
    wPlayback = nan(length(chans), length(indices));
    
    for index = indices
        
        if (strcmp(type,'m') || strcmp(type,'t')) && (index == 1)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_pos,...
                threshold,phase_at_0_pos,chans,desiredF(index),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_pos,fThresholdMin,fThresholdMax);
        elseif (strcmp(type,'m') || strcmp(type,'t')) && (index == 2)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_neg,...
                threshold,phase_at_0_neg,chans,desiredF(2),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_neg,fThresholdMin,fThresholdMax);
        elseif (strcmp(type,'s') && index ==1) || (strcmp(type,'t') && index == 3)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square,...
                threshold,phase_at_0,chans,desiredF,markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f,fThresholdMin,fThresholdMax);
        end
        
        count = 1;
        for i = chans
            mags = 1e6*dataForPPanalysis{i}{index}{1};
            label= dataForPPanalysis{i}{index}{4};
            keeps = dataForPPanalysis{i}{index}{5};
            maxLabel = max(unique(label));
            
            baseVal = avgFn(mags(label ==0 & keeps));
            difference = 100*(avgFn(mags(label ==3 & keeps)) - baseVal)/baseVal;
            percentInd = 100*(mags(label ==maxLabel & keeps) - baseVal)/baseVal;
            if baseVal > epThresholdMax
                wPlayback(count,index) = difference;
                wTotalPlayback(subjectNum,count,index) = difference;
                phaseTotal(subjectNum,count,index) = peakPhase(count);
            else
                wPlayBack(count,index) = nan;
                wTotalPlayback(subjectNum,count,index) = nan;
                phaseTotal(subjectNum,count,index) = nan;
            end
            count  = count + 1;
            
        end
        
        hPlayback(countScatter) =  scatter(peakPhase,wPlayback(:,index),markerSize,plotColor(subjectNum,:),'filled');
        
    end
    
    countScatter = countScatter + 1;
    
end
xlim([0 360])
xticks([0 45 90 135 180 225 270 315 360])
hline(0,'k')
legend(hPlayback,{'Subject 7',...
    'Subject 7 Playback'})
title('Phase of delivery and CEP modulation')
xlabel('Phase of delivery (degrees)')
ylabel({'Percent change in EP size','from baseline to >5 conditioning stimuli'})
set(gca,'fontsize',18)

%% null fit for subject 7
figure
clearvars hNull
hold on
countScatter = 1;

for sid = SIDS(end-1)
    sid = sid{:};
    subjid = sid;
    info = M(sid);
    type = info{1};
    subjectNum = info{3};
    desiredF = info{2};
    stims = info{4};
    bads = info{5};
    goodEPs = info{6};
    betaChan = info{7};
    chans = [1:64];
    badsTotal = [stims bads];
    chans(ismember(chans, badsTotal) | ~ismember(chans,goodEPs)) = [];
    Montage.MontageTokenized = {'Grid(1:64)'};
    
    load(fullfile(folderEP,strcat(subjid,['epSTATS-PP-sig' modifierEP '.mat'])))
    load(fullfile(folderPhase,[sid '_phaseDelivery_allChans' modifierPhase '.mat']));
    
    fprintf(['running for subject ' sid '\n']);
    
    %%
    indices = 3;
    
    for index = 1:2
        
        if (strcmp(type,'m') || strcmp(type,'t')) && (index == 1)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_pos,...
                threshold,phase_at_0_pos,chans,desiredF(index),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_pos,fThresholdMin,fThresholdMax);
        elseif (strcmp(type,'m') || strcmp(type,'t')) && (index == 2)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_neg,...
                threshold,phase_at_0_neg,chans,desiredF(2),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_neg,fThresholdMin,fThresholdMax);
        end
        wNull = nan(length(chans), indices);
        count = 1;
        for i = chans
            mags = 1e6*dataForPPanalysis{i}{index}{1};
            label= dataForPPanalysis{i}{index}{4};
            keeps = dataForPPanalysis{i}{index}{5};
            maxLabel = max(unique(label));
            
            baseVal = avgFn(mags(label ==0 & keeps));
            difference = 100*(avgFn(mags(label ==maxLabel & keeps)) - baseVal)/baseVal;
            if baseVal > epThresholdMax
                wNull(count,index) = difference;
                wTotalNull(subjectNum,count,index) = difference;
                phaseTotal(subjectNum,count,index) = peakPhase(count);
            else
                wNull(count,index) = nan;
                wTotalNull(subjectNum,count,index) = nan;
                phaseTotal(subjectNum,count,index) = nan;
            end
            count = count +1;
        end
        
        hNull(countScatter) =  scatter(peakPhase,wNull(:,index),markerSize,plotColor(subjectNum,:),'filled');
        
    end
    countScatter = countScatter + 1;
    
    index = 3;
    wNull = nan(length(chans), length(indices));
    count = 1;
    for i = chans
        mags = 1e6*dataForPPanalysis{i}{index}{1};
        label= dataForPPanalysis{i}{index}{4};
        keeps = dataForPPanalysis{i}{index}{5};
        baseVal = avgFn(mags(label ==0 & keeps));
        difference = 100*(avgFn(mags(label ==1 & keeps)) - baseVal)/baseVal;
        if baseVal > epThresholdMax
            wNull(count,index) = difference;
            wTotalNull(subjectNum,count,index) = difference;
            phaseTotal(subjectNum,count,index) = peakPhase(count);
        else
            wNull(count,index) = nan;
            wTotalNull(subjectNum,count,index) = nan;
            phaseTotal(subjectNum,count,index) = nan;
        end
        count = count +1;
    end
    
    hNull(countScatter) =  scatter(peakPhase,wNull(:,index),125,plotColor(subjectNum+2,:),'d','filled');
    
end
xlim([0 360])
ylim([-10 40])
hline(0,'k')
xticks([0 45 90 135 180 225 270 315 360])
legend(hNull,{'Subject 7',...
    'Subject 7 Null Control'})
title('Phase of delivery and CEP modulation')
xlabel('Phase of delivery (degrees)')
ylabel([{'Percent change in EP size from baseline'}])
set(gca,'fontsize',18)


%% all subjects fit all non-bad channels - see phases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
h = [];
for sid = SIDS
    
    sid = sid{:};
    subjid = sid;
    info = M(sid);
    type = info{1};
    subjectNum = info{3};
    desiredF = info{2};
    stims = info{4};
    bads = info{5};
    goodEPs = info{6};
    betaChan = info{7};
    chans = [1:64];
    badsTotal = [stims bads];
    chans(ismember(chans, badsTotal)) = [];
    Montage.MontageTokenized = {'Grid(1:64)'};
    
    
    load(fullfile(folderEP,strcat(subjid,['epSTATS-PP-sig' modifierEP '.mat'])))
    load(fullfile(folderPhase,[sid '_phaseDelivery_allChans' modifierPhase '.mat']));
    
    fprintf(['running for subject ' sid '\n']);
    
    %%
    
    if strcmp(type,'m')
        indices = [1,2];
    elseif strcmp(type,'s')
        indices = 1;
    elseif strcmp(type,'t')
        indices = [1,2,4];
    end
    
    for index = indices
        
        if (strcmp(type,'m') || strcmp(type,'t')) && (index == 1)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_pos,...
                threshold,phase_at_0_pos,chans,desiredF(index),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_pos,fThresholdMin,fThresholdMax);
            fDeliver = f_pos;
        elseif (strcmp(type,'m') || strcmp(type,'t')) && (index == 2)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_neg,...
                threshold,phase_at_0_neg,chans,desiredF(2),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_neg,fThresholdMin,fThresholdMax);
            fDeliver = f_neg;
            
        elseif (strcmp(type,'s') && index ==1) || (strcmp(type,'t') && index == 3)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square,...
                threshold,phase_at_0,chans,desiredF,markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f,fThresholdMin,fThresholdMax);
            fDeliver = f;
            
        end
        
        w = nan(length(chans), length(indices));
        count = 1;
        for i = chans
            w(count,index) = nanmean(fDeliver(:,i));
            count = count +1;
        end
        h(subjectNum) =  scatter(peakPhase,w(:,index),markerSize,plotColor(subjectNum,:),'filled');
        
    end
    
    
end
xlim([0 360])
xticks([0 45 90 135 180 225 270 315 360])
legend([h],{'Subject 1',...
    'Subject 2',...
    'Subject 3',...
    'Subject 4',...
    'Subject 5',...
    'Subject 6',...
    'Subject 7',...
    'Subject 7 Playback'})
title('Phase of delivery and frequency of fit beta')
xlabel('Phase of delivery (degrees)')
ylabel('Mean frequency of delivery')
set(gca,'fontsize',18)



%% subject 7 - with/without playback
% fit all non-bad channels - see phases playback only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
hPlayback = [];
countScatter = 1;
for sid = SIDS(end-1:end)
    
    sid = sid{:};
    subjid = sid;
    info = M(sid);
    type = info{1};
    subjectNum = info{3};
    desiredF = info{2};
    stims = info{4};
    bads = info{5};
    goodEPs = info{6};
    betaChan = info{7};
    chans = [1:64];
    badsTotal = [stims bads];
    chans(ismember(chans, badsTotal)) = [];
    Montage.MontageTokenized = {'Grid(1:64)'};
    
    
    load(fullfile(folderEP,strcat(subjid,['epSTATS-PP-sig' modifierEP '.mat'])))
    load(fullfile(folderPhase,[sid '_phaseDelivery_allChans' modifierPhase '.mat']));
    
    fprintf(['running for subject ' sid '\n']);
    
    %%
    
    if strcmp(type,'m')
        indices = [1,2];
    elseif strcmp(type,'s')
        indices = 1;
    elseif strcmp(type,'t')
        indices = [1,2,4];
    end
    
    for index = indices
        
        if (strcmp(type,'m') || strcmp(type,'t')) && (index == 1)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_pos,...
                threshold,phase_at_0_pos_acaus,chans,desiredF(index),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_pos,fThresholdMin,fThresholdMax);
            fDeliver = f_pos;
        elseif (strcmp(type,'m') || strcmp(type,'t')) && (index == 2)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square_neg,...
                threshold,phase_at_0_neg_acaus,chans,desiredF(2),markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f_neg,fThresholdMin,fThresholdMax);
            fDeliver = f_neg;
            
        elseif (strcmp(type,'s') && index ==1) || (strcmp(type,'t') && index == 3)
            [peakPhase,peakStd,peakLength,circularTest,markerSize] =  phase_delivery_accuracy_forPP(r_square,...
                threshold,phase_at_0,chans,desiredF,markerMin,markerMax,minData,maxData,markerToUse,testStatistic,f,fThresholdMin,fThresholdMax);
            fDeliver = f;
            
        end
        
        w = nan(length(chans), length(indices));
        count = 1;
        for i = chans
            w(count,index) = nanmean(fDeliver(:,i));
            count = count +1;
        end
        
        hPlayback(countScatter) =  scatter(peakPhase,w(:,index),markerSize,plotColor(subjectNum,:),'filled');
        
    end
    countScatter = countScatter + 1;
    
    
end
xlim([0 360])
xticks([0 45 90 135 180 225 270 315 360])
legend(hPlayback,{'Subject 7',...
    'Subject 7 Playback'})
title('Phase of delivery and frequency of fit beta')
xlabel('Phase of delivery (degrees)')
ylabel('Mean frequency of delivery')
set(gca,'fontsize',18)

