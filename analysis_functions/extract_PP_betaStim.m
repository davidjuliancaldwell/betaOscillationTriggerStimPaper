function [signalPP,pkLocs,trLocs] =  extract_PP_betaStim(signal,t,tBegin,tEnd,smooth,ppOpt,minPeakDist,framelen)
% extract_PP_betaStim
% extract peak to peak values in a signal
% this works on single trials
% time x trials
% ppOpt: 'abs' (default), 'falling', or 'rising' — passed to peak_to_peak_beta_stim
% minPeakDist: minimum distance between peaks in samples (default 15)
% framelen: Savitzky-Golay smoothing frame length (default 171; use 91 for
%           narrower extraction windows like 702d24)

% David.J.Caldwell
% 9.7.2018

if nargin < 6 || isempty(ppOpt)
    ppOpt = 'abs';
end
if nargin < 7 || isempty(minPeakDist)
    minPeakDist = 15;
end
if nargin < 8 || isempty(framelen)
    framelen = 171;
end

numTrials = size(signal,2);
signalPP = zeros(1,numTrials);
pkLocs = zeros(1,numTrials);
trLocs = zeros(1,numTrials);

order = 3;

plotIt = 0;

tTemp = t(t>tBegin & t<tEnd);
for i = 1:size(signal,2)

    tempSignalExtract = squeeze(signal(t>tBegin & t<tEnd,i));
    if smooth
        tempSignalExtract = sgolayfilt_complete(tempSignalExtract,order,framelen);
    end

    [amp,pk_loc,tr_loc]=peak_to_peak_beta_stim(tempSignalExtract,ppOpt,minPeakDist);
    
    if isempty(amp)
        amp = nan;
        pk_loc = nan;
        tr_loc = nan;
    end
    
    signalPP(i) = amp;
    pkLocs(i) = pk_loc;
    trLocs(i) = tr_loc;
    
    if plotIt && mod(i,10) == 0
        figure
        plot(tTemp,tempSignalExtract)
        vline(tTemp(pk_loc),'r')
        vline(tTemp(tr_loc),'b')
    end
    
end

end