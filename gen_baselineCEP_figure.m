% Regenerate smallMultiples baseline CEP figure from a saved .mat.
% Use case: re-do figures after fixing dependency issues without
% rerunning the expensive ECoG re-ref pass.
%
% Usage: set `sid` before sourcing, or call after setting idx in a loop.

setup_environment;
Z_Constants;

valueSet = {{'s',180,1,[54 62],[1 49 58 59],[44 45 46 52 53 55 60 61 63],53},...
    {'m',[0 180],2,[55 56],[1 2 3 31 57],[47 48 64],64},...
    {'s',180,3,[11 12],[57],[4 5 10 13],4},...
    {'s',270,4,[59 60],[1 9 10 35 43],[50 51 52 53 58],51},...
    {'m',[90,270],5,[13 14],[23 27 28 29 30 32 44 52 60],[5],5},...
    {'t',[270,90,12345,12345],6,[56 64],[57:63],[47 48 54 55],55},...
    {'m',[90,270],7,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31},...
    {'m',[90,270],8,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31}};
M = containers.Map(SIDS, valueSet, 'UniformValues', false);

if ~exist('sid', 'var')
    error('Set sid before running this script (e.g. sid = ''c91479'')');
end

info = M(sid);
stimChans = info{4};
betaChan = info{7};

matFile = fullfile(folderEP, [sid '_baselineCCEPs.mat']);
if ~exist(matFile, 'file')
    error('Missing %s. Run BETA_ExtractNeuralDataCEPscreen for %s first.', matFile, sid);
end
fprintf('Loading %s...\n', matFile);
S = load(matFile, 't', 'ECoGDataAverage');
t = S.t;
ECoGDataAverage = S.ECoGDataAverage;

smallMultiples(ECoGDataAverage, t, 'type1', stimChans, 'type2', betaChan, 'average', 1);
hFig = gcf;
set(hFig, 'Name', sprintf('%s — baseline CEP grid', sid));

outBase = fullfile(folderPlots, sprintf('baselineCEP_smallMultiples_%s', sid));
set(hFig, 'PaperPositionMode', 'auto');
print(hFig, [outBase '.png'], '-dpng', '-r300');
print(hFig, [outBase '.eps'], '-depsc', '-r600');
fprintf('Saved: %s.{png,eps}\n', outBase);
close(hFig);
