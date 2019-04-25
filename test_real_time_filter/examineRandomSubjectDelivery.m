%% try 702d24, 0b5a2e
%sid = input('what is the subject ID? ','s');
close all;clear all;clc

Z_Constants
saveIt = 0;

sid = '0b5a2e';
% c19479,7dbdec doesnt have continuous raw channel
switch sid
    case 'd5cd55'
        betaChan = 53;
        load(fullfile(folderECoGData,'d5cd55_forBetaPhase'))
        subject_num = '1';
        
    case '702d24'
        load(fullfile(folderECoGData,'702d24_forBetaPhase'))
        betaChan = 5;
        subject_num = '5';
        
    case 'c91479'
        load(fullfile(folderECoGData,'c91479_forBetaPhase'))
        betaChan = 64;
        subject_num = '2';
        
    case '0b5a2e'
        load(fullfile(folderECoGData,'0b5a2e_ECoG'))
        betaChan = 31;
        subject_num = '7';
        
    case '0b5a2ePlayback'
        betaChan = 31;
        load(fullfile(folderECoGData,'0b5a2e_Playback_ECoG'))
        
        subject_num = '7 playback';
        
    case '7dbdec'
        load(fullfile(folderECoGData,'7dbdec_forBetaPhase'))
        subject_num = '3';
        
        betaChan = 4;
    case 'ecb43e'
        load(fullfile(folderECoGData,'ecb43e_forBetaPhase'))
        subject_num = '6';
        
        betaChan = 55;
    case '9ab7ab'
        load(fullfile(folderECoGData,'9ab7ab_forBetaPhase'))
        betaChan = 51;
        subject_num = '4';
end

fprintf('loading in ecog data for %s:\n',sid);
fprintf('channel %d:\n',betaChan);
tic;
grp = floor((betaChan-1)/16);
ev = sprintf('ECO%d', grp+1);
achan = betaChan - grp*16;
formatSpec = 'ECO%d.data(:,%d)';
chanBlockInt = sprintf(formatSpec,grp+1,achan);

%         [eco, efs] = tdt_loadStream(tp, block, ev, achan);
raw_sig = 4.*eval(chanBlockInt);

filt_sig = 4.*Wave.data(:,3);
fs1 = ECO1.info.SamplingRateHz;
fs2 = Wave.info.SamplingRateHz;
stim_times = SMon.data(:,2);

% convert sampling rates
fac = fs2/fs1;

if (length(filt_sig)/2 ~= length(raw_sig)) & ~strcmp(sid,'7dbdec')
    filt_sig_decimate = [0; decimate(filt_sig,fac)]; % decimate because it's stored at double the rate of Eco
else
    filt_sig_decimate = decimate(filt_sig,fac);
end
%
t1 = 1e3*[0:length(filt_sig_decimate)-1]/fs1;

fig1 = figure;
plot(t1,raw_sig,'linewidth',2)
hold on
plot(t1,filt_sig_decimate,'linewidth',2)
timeStamps = find(stim_times>0);
timeStamps = 1e3*((timeStamps/2)/fs1);
vline([timeStamps],'k:');

legend({'Raw Signal','Filtered Signal','Stimulation Trigger'})
xlabel('time (ms)')
ylabel('amplitude')
set(gca,'fontsize', 14)
title('Operation of Real Time Filtering with Stimulation Blanking')
fig1.Units = "inches";
fig1.Position = [1 1 4 8];

clearvars -except saveIt raw_sig filt_sig_decimate stimTimes fac fs1 fs2 filt_sig t1 timeStamps sid betaChan subject_num folderData folderECoGData folderTiming folderPhase folderEP folderCoords folderPlots folderTestFilter

% here come the burst tables
if strcmp(sid,'0b5a2ePlayback')
    load(fullfile(folderTiming, ['0b5a2e_tables.mat']), 'bursts', 'fs', 'stims');
    delay = 577869;
else
    load(fullfile(folderTiming, [sid '_tables.mat']), 'bursts', 'fs', 'stims');
end

% get epoch's along burst
% adjust burst table to have 9th element which says which burst you are in
% for the conditioning bursts

% start off with nan, replace conditioning bursts with burst ID number
stims(9,:) = nan;

conditioning_epoched_filt = {};
conditioning_epoched_raw = {};

% presamp originally set to 0.020
preSamp = round(0.050 * fs1);
postSamp = round(0.05 * fs1);

utype = unique(bursts(5,:)); % unique stim types

% change sample id
stim_table_decimated = stims;
stim_table_decimated(2,:) = round(stims(2,:)/2);
t = (-preSamp:postSamp) / fs1;

if strcmp(sid,'d5cd55')
    bursts = bursts(:,(bursts(3,:)>4.5e6));
end


for ind = 1:size(bursts,2)
    % start index of bursts
    start = bursts(2,ind);
    % stop index of bursts
    stop = bursts(3,ind);
    
    % add in which burtst each stim pulse was in
    stims(9,(find((start<stims(2,:)) & (stims(2,:)<stop)))) = ind;
    cts = (stims(3,:)==1 & stims(9,:)==ind); % Stims that are conditioning stimuli
    
    conditioning_epoched_filt{ind} = squeeze(getEpochSignal(filt_sig_decimate, stim_table_decimated(2, cts)-preSamp, stim_table_decimated(2, cts)+postSamp+1)); %getting segments of beta signals
    conditioning_epoched_raw{ind} = squeeze(getEpochSignal(raw_sig, stim_table_decimated(2, cts)-preSamp, stim_table_decimated(2, cts)+postSamp+1)); %getting segments of raw ECoG signals
    
end

% mark stims that are null and plot them for 0b5a2e

if strcmp(sid,'0b5a2e')
    typeInt = 2;
    typei = 3; % this determines probe pulse to highlight
    
    % this is to modify burst table
    bursts_dec = bursts;
    bursts_dec(2,:) = round(bursts(2,:)/2);
    bursts_dec(3,:) = round(bursts(3,:)/2);
    
    bursts_ind_int = bursts(5,:) == typeInt;
    bursts_dec_sub = bursts_dec(:,bursts_ind_int);
    
    
    % highlight probe pulses
    pts = stims(3,:) == 0;
    pstims = stim_table_decimated(:,pts);
    types = unique(bursts(5,pstims(4,:)));
    
    probes = pstims(5,:) < 0.250*fs & bursts(5,pstims(4,:))==types(typei);
    probes = pstims(5,:) < 0.250*fs;
    probe_times = pstims(:,probes);
    
    % conditioning
    cond_pts = stims(3,:) ==1;
    cond_pstims = stim_table_decimated(:,cond_pts);
    
    
    fig = figure;
    plot(t1,1e6*filt_sig_decimate,'linewidth',2)
    % vline([1e3*probe_times(2,:)/fs1],'k:');
    vline([1e3*pstims(2,:)/fs1],'k:')
    vline([1e3*cond_pstims(2,:)/fs1],'r:');
    %legend({'Filtered Signal','Null Probes','Cond Probes'})
    for i = 1:size(bursts_dec_sub,2)  
        highlight(gca, [1e3*bursts_dec_sub(2,i)/fs1 1e3*bursts_dec_sub(3,i)/fs1], [], [.5 .5 .5]) %this is the part that plots that stim window     
    end
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
    xlim([4.79e5 4.83e5])
    set(gca,'fontsize', 14)
    title('Operation of Real Time Filtering with Stimulation Blanking')
    fig.Units = "inches";
    fig.Position = [0 0 8 4];
    %
end

% look at bursts for stavros that are greate than five.
%%

% bursts == 3 for ecb43e, otherwise 0,1, for 0b5a2e
[indices] = find((bursts(4,:) >= 5) & bursts(5,:) == 0);

indices_rand = datasample(indices,4,'Replace',false);
%
fig1 = figure;
fig1.Units = 'inches';
fig1.Position =   [0 0 4 4];

for i = 1:4
    subplot(2,2,i)
    hold on
    axis tight
    
    plot(1e3*t,1e6*conditioning_epoched_filt{indices_rand(i)})
    plot(1e3*t,1e6*mean(conditioning_epoched_filt{indices_rand(i)},2),'k','linewidth',2)
    vline(0)
    set(gca,'fontsize',10)
end
xlabel('Time (ms)')
ylabel(['Amplitude (\mu V)'])
%subtitle(['Filtered beta signal during beta burst for subject ' subject_num]);
%title('Filtered Signal')
if saveIt
    SaveFig(folderPlots, sprintf(['betaBurst-subj-%s-v3'], subject_num), 'eps', '-r600');
    SaveFig(folderPlots, sprintf(['betaBurst-subj-%s-v3'], subject_num), 'png', '-r600');
    SaveFig(folderPlots, sprintf(['betaBurst-subj-%s-v3'], subject_num), 'svg', '-r600');
end

for i = 1:4
    subplot(2,2,i)
    hold on
    axis tight
    
    plot(1e3*t,conditioning_epoched_filt{indices_rand(i)})
    plot(1e3*t,mean(conditioning_epoched_filt{indices_rand(i)},2),'k','linewidth',2)
    vline(0)
    set(gca,'fontsize',12)
end
%% visualize

% interest
for i = 1:4
    figure
    subplot(2,1,1)
    hold on
    plot(1e3*t,conditioning_epoched_filt{indices_rand(i)})
    plot(1e3*t,mean(conditioning_epoched_filt{indices_rand(i)},2),'k','linewidth',2)
    vline(0)
    xlabel('time (ms)')
    ylabel('amplitude')
    title('Filtered Signal')
    
    subplot(2,1,2)
    hold on
    plot(1e3*t,conditioning_epoched_raw{indices_rand(i)})
    plot(1e3*t,mean(conditioning_epoched_raw{indices_rand(i)},2),'k','linewidth',2)
    ylim([-1e-4 1e-4])
    vline(0)
    xlabel('time (ms)')
    ylabel('amplitude')
    title('Raw Signal')
end


