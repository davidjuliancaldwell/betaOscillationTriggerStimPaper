% Inspect variable structure in phase data files
phaseDir = '/Users/davidcaldwell/code/betaOscillationTriggerStimPaper/data/phase_data/';

fprintf('\n=== Multi-phase subject: c91479 ===\n');
S = load(fullfile(phaseDir, 'c91479_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat'));
disp(fieldnames(S));
fns = fieldnames(S);
for i = 1:length(fns)
    v = S.(fns{i});
    fprintf('  %s: class=%s, size=[%s]\n', fns{i}, class(v), num2str(size(v)));
end

fprintf('\n=== Single-phase subject: d5cd55 ===\n');
S2 = load(fullfile(phaseDir, 'd5cd55_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat'));
disp(fieldnames(S2));
fns2 = fieldnames(S2);
for i = 1:length(fns2)
    v = S2.(fns2{i});
    fprintf('  %s: class=%s, size=[%s]\n', fns2{i}, class(v), num2str(size(v)));
end

fprintf('\n=== Multi-phase subject: ecb43e ===\n');
S3 = load(fullfile(phaseDir, 'ecb43e_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat'));
disp(fieldnames(S3));
fns3 = fieldnames(S3);
for i = 1:length(fns3)
    v = S3.(fns3{i});
    fprintf('  %s: class=%s, size=[%s]\n', fns3{i}, class(v), num2str(size(v)));
end

% Check range of phase values to see if radians or degrees
fprintf('\n=== Phase range check (c91479, chan 39) ===\n');
fprintf('  phase_at_0_pos(39,:) min=%.3f max=%.3f\n', min(S.phase_at_0_pos(39,:)), max(S.phase_at_0_pos(39,:)));
fprintf('  phase_at_0_neg(39,:) min=%.3f max=%.3f\n', min(S.phase_at_0_neg(39,:)), max(S.phase_at_0_neg(39,:)));

fprintf('\n=== Phase range check (d5cd55, chan 44) ===\n');
fprintf('  phase_at_0(44,:) min=%.3f max=%.3f\n', min(S2.phase_at_0(44,:)), max(S2.phase_at_0(44,:)));

% Check r_square and frequency structure
fprintf('\n=== r_square and frequency check (c91479) ===\n');
fprintf('  r_square_pos(39,:) first 5: '); fprintf('%.3f ', S.r_square_pos(39,1:5)); fprintf('\n');
fprintf('  frequency_pos(39,:) first 5: '); fprintf('%.3f ', S.frequency_pos(39,1:5)); fprintf('\n');
