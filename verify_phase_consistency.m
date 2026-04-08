% verify_phase_consistency.m
% Check that phase fit data in .mat files is consistent with expected phase labels.
% For each subject, load phase data, filter to good fits (r_square>0.7, freq 12.01-19.99),
% and report circular mean phase, mean r_square, mean frequency per EP channel.
% Uses circ_mean from CircStat toolbox (returns value in [-pi, pi]).

phaseDir = fullfile(pwd, 'data', 'phase_data');

% Add CircStat toolbox to path
addpath(genpath('/Users/davidcaldwell/code/matlab_utilities/circular_statistics/CircStat2012a'));

% Helper: convert circ_mean output (radians, [-pi,pi]) to degrees [0, 360)
% circ_mean returns in [-pi, pi], rad2deg gives [-180,180], mod wraps to [0,360)
circ_mean_deg = @(phases_rad) mod(rad2deg(circ_mean(phases_rad)), 360);

%% EP channels per subject (from valueSet in multipleSubj_GLMM_script_PP.m)
epChans.d5cd55 = [44 45 46 52 53 55 60 61 63];
epChans.c91479 = [47 48 64];
epChans.x7dbdec = [4 5 10 13];
epChans.x9ab7ab = [50 51 52 53 58];
epChans.x702d24 = [5];
epChans.ecb43e = [47 48 54 55];
epChans.x0b5a2e = [14 15 16 20 21 23 31 32 40];
epChans.x0b5a2ePlayBack = [14 15 16 20 21 23 31 32 40];

%% Beta reference channel per subject (from valueSet field 7)
betaRef.d5cd55 = 53;
betaRef.c91479 = 64;
betaRef.x7dbdec = 4;
betaRef.x9ab7ab = 51;
betaRef.x702d24 = 5;
betaRef.ecb43e = 55;
betaRef.x0b5a2e = 31;
betaRef.x0b5a2ePlayBack = 31;

%% Expected phase targets (degrees)
expected.c91479 = struct('pos_target', 0, 'neg_target', 180, 'desiredF', '[0, 180]');
expected.x702d24 = struct('pos_target', 90, 'neg_target', 270, 'desiredF', '[90, 270]');
expected.ecb43e = struct('pos_target', 270, 'neg_target', 90, 'desiredF', '[270, 90, rand, rand]');
expected.x0b5a2e = struct('pos_target', 90, 'neg_target', 270, 'desiredF', '[90, 270]');
expected.x0b5a2ePlayBack = struct('pos_target', 90, 'neg_target', 270, 'desiredF', '[90, 270]');
expected.d5cd55 = struct('target', 180, 'desiredF', '180');
expected.x7dbdec = struct('target', 180, 'desiredF', '180');
expected.x9ab7ab = struct('target', 270, 'desiredF', '270');

%% ========================================
%  MULTI-PHASE SUBJECTS
%  ========================================
multiSubjects = {'c91479', '702d24', 'ecb43e', '0b5a2e', '0b5a2ePlayBack'};
multiFieldNames = {'c91479', 'x702d24', 'ecb43e', 'x0b5a2e', 'x0b5a2ePlayBack'};

for si = 1:length(multiSubjects)
    sid = multiSubjects{si};
    fn = multiFieldNames{si};

    fname = fullfile(phaseDir, [sid '_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat']);
    fprintf('\n=========================================================\n');
    fprintf('MULTI-PHASE SUBJECT: %s\n', sid);
    fprintf('Expected: desiredF = %s, pos -> %d deg, neg -> %d deg\n', ...
        expected.(fn).desiredF, expected.(fn).pos_target, expected.(fn).neg_target);
    fprintf('=========================================================\n');

    S = load(fname);
    nChan = size(S.phase_at_0_pos, 2);
    channels = epChans.(fn);

    fprintf('\n--- POS condition (target: %d deg) ---\n', expected.(fn).pos_target);
    fprintf('%-6s %-12s %-12s %-12s %-10s %-10s %-10s\n', ...
        'Chan', 'CircMean(d)', 'MeanRsq', 'MeanFreq', 'nTrials', 'nGood', 'pctGood');

    for ch = channels
        if ch > nChan
            fprintf('%-6d  *** Channel %d exceeds data size (%d) ***\n', ch, ch, nChan);
            continue;
        end

        phases = S.phase_at_0_pos(:, ch);
        rsq = S.r_square_pos(:, ch);
        freq = S.f_pos(:, ch);

        goodMask = (rsq > 0.7) & (freq > 12.01) & (freq < 19.99);
        nTotal = length(phases);
        nGood = sum(goodMask);

        if nGood > 0
            cm = circ_mean_deg(phases(goodMask));
            mr = mean(rsq(goodMask));
            mf = mean(freq(goodMask));
        else
            cm = NaN; mr = NaN; mf = NaN;
        end

        fprintf('%-6d %-12.1f %-12.3f %-12.2f %-10d %-10d %-10.1f%%\n', ...
            ch, cm, mr, mf, nTotal, nGood, 100*nGood/nTotal);
    end

    fprintf('\n--- NEG condition (target: %d deg) ---\n', expected.(fn).neg_target);
    fprintf('%-6s %-12s %-12s %-12s %-10s %-10s %-10s\n', ...
        'Chan', 'CircMean(d)', 'MeanRsq', 'MeanFreq', 'nTrials', 'nGood', 'pctGood');

    for ch = channels
        if ch > nChan
            fprintf('%-6d  *** Channel %d exceeds data size (%d) ***\n', ch, ch, nChan);
            continue;
        end

        phases = S.phase_at_0_neg(:, ch);
        rsq = S.r_square_neg(:, ch);
        freq = S.f_neg(:, ch);

        goodMask = (rsq > 0.7) & (freq > 12.01) & (freq < 19.99);
        nTotal = length(phases);
        nGood = sum(goodMask);

        if nGood > 0
            cm = circ_mean_deg(phases(goodMask));
            mr = mean(rsq(goodMask));
            mf = mean(freq(goodMask));
        else
            cm = NaN; mr = NaN; mf = NaN;
        end

        fprintf('%-6d %-12.1f %-12.3f %-12.2f %-10d %-10d %-10.1f%%\n', ...
            ch, cm, mr, mf, nTotal, nGood, 100*nGood/nTotal);
    end

    % Consistency check
    fprintf('\n--- CONSISTENCY CHECK ---\n');
    pos_target = expected.(fn).pos_target;
    neg_target = expected.(fn).neg_target;

    for ch = channels
        if ch > nChan, continue; end

        phases_p = S.phase_at_0_pos(:, ch);
        rsq_p = S.r_square_pos(:, ch);
        freq_p = S.f_pos(:, ch);
        gm_p = (rsq_p > 0.7) & (freq_p > 12.01) & (freq_p < 19.99);

        phases_n = S.phase_at_0_neg(:, ch);
        rsq_n = S.r_square_neg(:, ch);
        freq_n = S.f_neg(:, ch);
        gm_n = (rsq_n > 0.7) & (freq_n > 12.01) & (freq_n < 19.99);

        if sum(gm_p) > 0 && sum(gm_n) > 0
            cm_p = circ_mean_deg(phases_p(gm_p));
            cm_n = circ_mean_deg(phases_n(gm_n));

            diff_p = abs(mod(cm_p - pos_target + 180, 360) - 180);
            diff_n = abs(mod(cm_n - neg_target + 180, 360) - 180);

            status_p = 'OK'; if diff_p > 90, status_p = '*** MISMATCH ***'; end
            status_n = 'OK'; if diff_n > 90, status_n = '*** MISMATCH ***'; end

            fprintf('  Chan %d: pos circMean=%.1f (target=%d, diff=%.1f) %s | neg circMean=%.1f (target=%d, diff=%.1f) %s\n', ...
                ch, cm_p, pos_target, diff_p, status_p, cm_n, neg_target, diff_n, status_n);
        else
            fprintf('  Chan %d: insufficient good fits (pos: %d, neg: %d)\n', ch, sum(gm_p), sum(gm_n));
        end
    end
end

%% ========================================
%  c91479: VERIFY 0 vs 180 TARGET ASSIGNMENT
%  ========================================
fprintf('\n=========================================================\n');
fprintf('c91479: VERIFY 0 vs 180 TARGET ASSIGNMENT\n');
fprintf('  desiredF(1)=0 assigned to pos (stims8==1)\n');
fprintf('  desiredF(2)=180 assigned to neg (stims8==0)\n');
fprintf('=========================================================\n');

S_c91 = load(fullfile(phaseDir, ...
    'c91479_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat'));
channels_c91 = epChans.c91479;
nChan_c91 = size(S_c91.phase_at_0_pos, 2);

fprintf('\n%-6s | %-14s %-14s | %-14s %-14s | %s\n', ...
    'Chan', 'pos circMean', 'closer to', 'neg circMean', 'closer to', 'Verdict');
fprintf('%s\n', repmat('-', 1, 90));

for ch = channels_c91
    if ch > nChan_c91
        fprintf('%-6d | exceeds data matrix (%d cols)\n', ch, nChan_c91);
        continue;
    end

    phases_p = S_c91.phase_at_0_pos(:, ch);
    rsq_p    = S_c91.r_square_pos(:, ch);
    freq_p   = S_c91.f_pos(:, ch);
    gm_p     = (rsq_p > 0.7) & (freq_p > 12.01) & (freq_p < 19.99);

    phases_n = S_c91.phase_at_0_neg(:, ch);
    rsq_n    = S_c91.r_square_neg(:, ch);
    freq_n   = S_c91.f_neg(:, ch);
    gm_n     = (rsq_n > 0.7) & (freq_n > 12.01) & (freq_n < 19.99);

    if sum(gm_p) == 0 || sum(gm_n) == 0
        fprintf('%-6d | insufficient good fits (pos:%d, neg:%d)\n', ch, sum(gm_p), sum(gm_n));
        continue;
    end

    cm_p = circ_mean_deg(phases_p(gm_p));
    cm_n = circ_mean_deg(phases_n(gm_n));

    dist_p_to_0   = abs(mod(cm_p - 0   + 180, 360) - 180);
    dist_p_to_180 = abs(mod(cm_p - 180 + 180, 360) - 180);
    dist_n_to_0   = abs(mod(cm_n - 0   + 180, 360) - 180);
    dist_n_to_180 = abs(mod(cm_n - 180 + 180, 360) - 180);

    if dist_p_to_0 < dist_p_to_180, closer_p = '0'; else, closer_p = '180'; end
    if dist_n_to_0 < dist_n_to_180, closer_n = '0'; else, closer_n = '180'; end

    if strcmp(closer_p, '0') && strcmp(closer_n, '180')
        verdict = 'CORRECT (pos~0, neg~180)';
    elseif strcmp(closer_p, '180') && strcmp(closer_n, '0')
        verdict = '*** SWAPPED (pos~180, neg~0) ***';
    else
        verdict = sprintf('MIXED (pos~%s, neg~%s)', closer_p, closer_n);
    end

    fprintf('%-6d | %-14.1f %-14s | %-14.1f %-14s | %s\n', ...
        ch, cm_p, closer_p, cm_n, closer_n, verdict);
end

%% ========================================
%  SINGLE-PHASE SUBJECTS: TARGET vs MEASURED
%  ========================================
% For each single-phase subject, verify the measured circular mean
% is closer to the stated target than to the opposite phase (target+180).

singleSubjects = {'d5cd55', '7dbdec', '9ab7ab'};
singleFieldNames = {'d5cd55', 'x7dbdec', 'x9ab7ab'};

for si = 1:length(singleSubjects)
    sid = singleSubjects{si};
    fn = singleFieldNames{si};
    target = expected.(fn).target;
    opposite = mod(target + 180, 360);

    fname = fullfile(phaseDir, [sid '_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat']);
    fprintf('\n=========================================================\n');
    fprintf('SINGLE-PHASE: %s  (target=%d, opposite=%d)\n', sid, target, opposite);
    fprintf('=========================================================\n');

    S = load(fname);
    nChan = size(S.phase_at_0, 2);
    channels = epChans.(fn);

    fprintf('\n%-6s | %-14s %-14s %-14s | %s\n', ...
        'Chan', 'circMean', 'distToTarget', 'distToOpposite', 'Verdict');
    fprintf('%s\n', repmat('-', 1, 75));

    for ch = channels
        if ch > nChan
            fprintf('%-6d | exceeds data matrix (%d cols)\n', ch, nChan);
            continue;
        end

        phases = S.phase_at_0(:, ch);
        rsq = S.r_square(:, ch);
        freq = S.f(:, ch);
        gm = (rsq > 0.7) & (freq > 12.01) & (freq < 19.99);

        if sum(gm) == 0
            fprintf('%-6d | no good fits (%d total trials)\n', ch, length(phases));
            continue;
        end

        cm = circ_mean_deg(phases(gm));
        dist_target   = abs(mod(cm - target   + 180, 360) - 180);
        dist_opposite = abs(mod(cm - opposite + 180, 360) - 180);

        if dist_target < dist_opposite
            verdict = sprintf('CORRECT (closer to %d)', target);
        else
            verdict = sprintf('*** CLOSER TO OPPOSITE %d ***', opposite);
        end

        fprintf('%-6d | %-14.1f %-14.1f %-14.1f | %s\n', ...
            ch, cm, dist_target, dist_opposite, verdict);
    end
end

%% ========================================
%  ADDITIONAL: ecb43e "random" condition
%  ========================================
fprintf('\n=========================================================\n');
fprintf('ECBAE: RANDOM condition (phase_at_0, non-targeted blocks)\n');
fprintf('=========================================================\n');
S = load(fullfile(phaseDir, 'ecb43e_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat'));
channels = epChans.ecb43e;
nChan = size(S.phase_at_0, 2);
fprintf('%-6s %-12s %-12s %-12s %-10s %-10s\n', ...
    'Chan', 'CircMean(d)', 'MeanRsq', 'MeanFreq', 'nTrials', 'nGood');
for ch = channels
    if ch > nChan, continue; end
    phases = S.phase_at_0(:, ch);
    rsq = S.r_square(:, ch);
    freq = S.f(:, ch);
    goodMask = (rsq > 0.7) & (freq > 12.01) & (freq < 19.99);
    nTotal = length(phases);
    nGood = sum(goodMask);
    if nGood > 0
        cm = circ_mean_deg(phases(goodMask));
        mr = mean(rsq(goodMask));
        mf = mean(freq(goodMask));
    else
        cm = NaN; mr = NaN; mf = NaN;
    end
    fprintf('%-6d %-12.1f %-12.3f %-12.2f %-10d %-10d\n', ...
        ch, cm, mr, mf, nTotal, nGood);
end

%% ========================================
%  BETA REFERENCE CHANNEL SUMMARY
%  ========================================
% For each subject, show the beta reference channel's measured phase
% vs target, and whether it maps to the correct phaseClass bin.
% Binning rule: circMean > 180 → 270, circMean <= 180 → 90.

fprintf('\n=========================================================\n');
fprintf('BETA REFERENCE CHANNEL SUMMARY\n');
fprintf('=========================================================\n');

allSubjects = {'d5cd55','c91479','7dbdec','9ab7ab','702d24','ecb43e','0b5a2e','0b5a2ePlayBack'};
allFieldNames = {'d5cd55','c91479','x7dbdec','x9ab7ab','x702d24','ecb43e','x0b5a2e','x0b5a2ePlayBack'};

fprintf('\n%-16s %-8s %-6s %-10s %-10s %-10s %-6s %s\n', ...
    'Subject', 'betaCh', 'Cond', 'Target', 'CircMean', 'Diff', 'Bin', 'Status');
fprintf('%s\n', repmat('-', 1, 85));

for si = 1:length(allSubjects)
    sid = allSubjects{si};
    fn = allFieldNames{si};
    bch = betaRef.(fn);

    fname = fullfile(phaseDir, [sid '_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat']);
    if ~exist(fname, 'file'), continue; end
    S = load(fname);

    if isfield(expected.(fn), 'pos_target')
        % Multi-phase subject
        nChan = size(S.phase_at_0_pos, 2);
        if bch > nChan
            fprintf('%-16s %-8d  *** channel exceeds data (%d) ***\n', sid, bch, nChan);
            continue;
        end

        % pos
        gm = (S.r_square_pos(:,bch) > 0.7) & (S.f_pos(:,bch) > 12.01) & (S.f_pos(:,bch) < 19.99);
        if sum(gm) > 0
            cm = circ_mean_deg(S.phase_at_0_pos(gm, bch));
            tgt = expected.(fn).pos_target;
            d = abs(mod(cm - tgt + 180, 360) - 180);
            bin = 270; if cm <= 180, bin = 90; end
            opp = mod(tgt + 180, 360);
            d_opp = abs(mod(cm - opp + 180, 360) - 180);
            if d < d_opp, st = 'OK'; else, st = '*** WRONG HALF ***'; end
            fprintf('%-16s %-8d %-6s %-10d %-10.1f %-10.1f %-6d %s\n', sid, bch, 'pos', tgt, cm, d, bin, st);
        end

        % neg
        gm = (S.r_square_neg(:,bch) > 0.7) & (S.f_neg(:,bch) > 12.01) & (S.f_neg(:,bch) < 19.99);
        if sum(gm) > 0
            cm = circ_mean_deg(S.phase_at_0_neg(gm, bch));
            tgt = expected.(fn).neg_target;
            d = abs(mod(cm - tgt + 180, 360) - 180);
            bin = 270; if cm <= 180, bin = 90; end
            opp = mod(tgt + 180, 360);
            d_opp = abs(mod(cm - opp + 180, 360) - 180);
            if d < d_opp, st = 'OK'; else, st = '*** WRONG HALF ***'; end
            fprintf('%-16s %-8d %-6s %-10d %-10.1f %-10.1f %-6d %s\n', sid, bch, 'neg', tgt, cm, d, bin, st);
        end
    else
        % Single-phase subject
        nChan = size(S.phase_at_0, 2);
        if bch > nChan
            fprintf('%-16s %-8d  *** channel exceeds data (%d) ***\n', sid, bch, nChan);
            continue;
        end
        gm = (S.r_square(:,bch) > 0.7) & (S.f(:,bch) > 12.01) & (S.f(:,bch) < 19.99);
        if sum(gm) > 0
            cm = circ_mean_deg(S.phase_at_0(gm, bch));
            tgt = expected.(fn).target;
            d = abs(mod(cm - tgt + 180, 360) - 180);
            bin = 270; if cm <= 180, bin = 90; end
            opp = mod(tgt + 180, 360);
            d_opp = abs(mod(cm - opp + 180, 360) - 180);
            if d < d_opp, st = 'OK'; else, st = '*** WRONG HALF ***'; end
            fprintf('%-16s %-8d %-6s %-10d %-10.1f %-10.1f %-6d %s\n', sid, bch, 'all', tgt, cm, d, bin, st);
        end
    end
end

fprintf('\n\n=== DONE ===\n');
