function [] = plot_brains_peak_func_avg(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,indices,saveFig,OUTPUT_DIR)
%% plot differences

epThresholdMin = 25;
epThresholdMax = 1500;
epThresholdAverage = 50;

cmap = cbrewer('seq','Purples',40);

w = nan(size(Grid, 1), 1);

for i = 1:64
    labelTotal = [];
    keepsTotal = [];
    magsTotal = [];
    if ~any(i==badsTotal) && any(i ==goodEPs)
        for index = indices
            mags = 1e6*dataForPPanalysis{i}{index}{1};
            mags(mags<epThresholdMin) = nan;
            mags(mags>epThresholdMax) = nan;
            label= dataForPPanalysis{i}{index}{4};
            keeps = dataForPPanalysis{i}{index}{5};
            
            labelTotal = [labelTotal label];
            keepsTotal = [keepsTotal keeps];
            magsTotal = [magsTotal mags];
        end
        maxLabel = max(unique(labelTotal));
        ppMax = nanmean(magsTotal(labelTotal ==maxLabel & keepsTotal));
        w(i) = ppMax;
        ppMax
    end
end

clims = [0 max(w)];

figure
set(gcf, 'Units', 'pixels', 'OuterPosition', [286.6000 108.6000 788.0000 645.6000]);

% plot beta channel overlaid
betaChanPlot = PlotBrainJustDots(subjid,{betaChan},[255, 153, 0]/255,true,800);
PlotDotsDirect(subjid, Grid, w, determineHemisphereOfCoverage(subjid), clims, 20, cmap, 1:size(Grid, 1), true);

% plot stimulation channels
stimulationPlot = PlotBrainJustDots(subjid,{stims(1),stims(2)},[0 0 0; 0 0 0],true);
leg = legend([stimulationPlot(1),stimulationPlot(2),betaChanPlot],...
    {['stimulation channel'],['stimulation channel'],...
    ['trigger channel = ' num2str(betaChan)]},'location','southwest');

colormap(cmap);
h = colorbar;
ylabel(h,'Volts (\muV)')
if subjectNum == 8
    title({['Subject 7 Playback'], ['Peak to peak evoked potential magnitude ']})
    
else
    title({['Subject ' num2str(subjectNum)], 'Peak to peak evoked potential magnitude '})
end
set(gca,'fontsize', 14)

if saveFig
    %     SaveFig(OUTPUT_DIR, sprintf(['EP-phase-%d-sid-%s-chan-%d'],typei,sid, chan,type,signalType), 'svg');
    SaveFig(OUTPUT_DIR, sprintf(['cortex-EP-phase-%d-sid-%s'],index,sid), 'png','-r300');
    close
end

%% plot percent differences
cmap = flipud(cbrewer('div','PiYG',40));

w = nan(size(Grid, 1), 1);

for i = 1:64
    labelTotal = [];
    keepsTotal = [];
    magsTotal = [];
    if ~any(i==badsTotal) && any(i ==goodEPs)
        for index = indices
            mags = 1e6*dataForPPanalysis{i}{index}{1};
            
            mags(mags<epThresholdMin) = nan;
            mags(mags>epThresholdMax) = nan;
            
            label= dataForPPanalysis{i}{index}{4};
            keeps = dataForPPanalysis{i}{index}{5};
            
            labelTotal = [labelTotal label];
            keepsTotal = [keepsTotal keeps];
            magsTotal = [magsTotal mags];
        end
        difference = 100*(nanmean(magsTotal(labelTotal ==3 & keepsTotal)) - nanmean(magsTotal(labelTotal ==0 & keepsTotal)))/nanmean(magsTotal(labelTotal ==0 & keepsTotal));
        
        if nanmean(magsTotal(labelTotal ==0 & keepsTotal)) > epThresholdAverage
            w(i) = difference;
        end
    end
end

clims = [-max(abs(min(w)),abs(max(w))) max(abs(min(w)),abs(max(w)))];

figure
set(gcf, 'Units', 'pixels', 'OuterPosition', [286.6000 108.6000 788.0000 645.6000]);

betaChanPlot = PlotBrainJustDots(subjid,{betaChan},[255, 153, 0]/255,true,800);

PlotDotsDirect(subjid, Grid, w, determineHemisphereOfCoverage(subjid), clims, 20, cmap, 1:size(Grid, 1), true);

% plot stimulation channels
stimulationPlot = PlotBrainJustDots(subjid,{stims(1),stims(2)},[0 0 0; 0 0 0],true);

leg = legend([stimulationPlot(1),stimulationPlot(2),betaChanPlot],...
    {['stimulation channel'],['stimulation channel'],...
    ['trigger channel = ' num2str(betaChan)]},'location','southwest');
colormap(cmap);
h = colorbar;
if subjectNum == 8
    title({['Subject 7 Playback'], ['Percent Baseline and Test Pulse (>5 stims) CEP difference']})
    
else
    title({['Subject ' num2str(subjectNum)], 'Percent Baseline and Test Pulse (>5 stims) CEP difference'})
end
ylabel(h,'Percent Difference')
set(gca,'fontsize', 14)

if saveFig
    SaveFig(OUTPUT_DIR, sprintf(['cortex-percentChange-EP-phase-%d-sid-%s'],index,sid), 'png','-r300');
    close
end
end