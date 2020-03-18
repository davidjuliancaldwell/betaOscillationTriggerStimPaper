%% plot peak to peak differences on individual cortical surfaces

setup_environment
%% parameters
% user defined
modifier = '-reref-50-new';
plotAvg = 1;
plotByAprioriPhase = 0;
saveFig = 0;
SIDS = {'d5cd55','c91479','7dbdec','9ab7ab','702d24','ecb43e','0b5a2e','0b5a2ePlayback'};


% valueSet = {{'s',180,1,[54 62],[1 49 58 59],[44 45 46 47 48 52 53 55 60 61 63],53},...
%     {'m',[0 180],2,[55 56],[1 2 3 31 57],[31 39 40 47 48 63 64],64},...
%     {'s',180,3,[11 12],[57],[4 5 10 13 18 19 20],4},...
%     {'s',270,4,[59 60],[1 9 10 35 43],[41 42 43 44 45 49 50 51 52 53 57 58 61 62],51},...
%     {'m',[90,270],5,[13 14],[23 27 28 29 30 32 44 52 60],[5],5},...
%     {'t',[90,180],6,[56 64],[57:64],[46 48 54 55 63],55},...
%     {'m',[90,270],7,[22 30],[24 25 29],[13 14 15 16 20 21 23 31 32 39 40],31},...
%     {'m',[90,270],8,[22 30],[24 25 29],[13 14 15 16 20 21 23 31 32 39 40],31}};

valueSet = {{'s',180,1,[54 62],[1 49 58 59],[44 45 46 52 53 55 60 61 63],53,2.5},...
    {'m',[0 180],2,[55 56],[1 2 3 31 57],[47 48 64],64,3},...
    {'s',180,3,[11 12],[57],[4 5 10 13],4,3.5},...
    {'s',270,4,[59 60],[1 9 10 35 43],[50 51 52 53 58],51,0.75},...
    {'m',[90,270],5,[13 14],[23 27 28 29 30 32 44 52 60],[5],5,0.75},...
    {'t',[270,90,12345,12345],6,[56 64],[57:63],[47 48 54 55 63],55,1.75}...
    {'m',[90,270],7,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31,1.75},...
    {'m',[90,270],8,[22 30],[24 25 29],[14 15 16 20 21 23 31 32 40],31,1.75}};

M = containers.Map(SIDS,valueSet,'UniformValues',false);

SIDS = {'c91479','9ab7ab'};

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
    chans(ismember(chans, badsTotal) | ~ismember(chans,goodEPs)) = [];
    Montage.MontageTokenized = {'Grid(1:64)'};
    
    load(strcat(subjid,['epSTATS-PP-sig' modifier '.mat']))
    
    if (strcmpi(sid,'0b5a2ePlayback'))
        load(fullfile(folderCoords,'0b5a2e','trodes.mat'));
        subjid = '0b5a2e';
    else
        load(fullfile(folderCoords,sid,'trodes.mat'))
    end
    
    %% plotting average deflection on all PP amplitudes regardless of phase
    
    if plotAvg
        if strcmp(type,'s')
            index = 1;
            plot_brains_peak_func_avg(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,index,saveFig,folderPlots)
            
        elseif strcmp(type,'m')
            
            index = 1:2;
            plot_brains_peak_func_avg(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,index,saveFig,folderPlots)
            
        elseif strcmp(type,'t')
            index = [1,2,4] ;
            plot_brains_peak_func_avg(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,index,saveFig,folderPlots)
            
        end
    end
    
    %% plotting average deflection on all PP amplitudes by apriori phase
    if plotByAprioriPhase
        if strcmp(type,'s')
            index = 1;
            plot_brains_peak_func(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,index,saveFig,folderPlots)
            
        elseif strcmp(type,'m')
            
            for index = 1:2
                plot_brains_peak_func(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,index,saveFig,folderPlots)
            end
            
        elseif strcmp(type,'t')
            for index = [1,2,4]
                plot_brains_peak_func_avg(dataForPPanalysis,subjid,sid,subjectNum,Grid,betaChan,stims,badsTotal,goodEPs,index,saveFig,folderPlots)
            end
            
        end
    end
end

