function tbt=makeTbtForDistractorLED(alignmentDir, alignment, onlyLongITI, jitterInSec)

if ~isfield(alignment,'isGrooming')
    alignment.isGrooming=zeros(size(alignment.cue));
end

if onlyLongITI==1
    % only take LED distractors occuring after an ITI similar to real cue
    % ITI
    cueThresh=0.5;
    distractThresh=0.5;
    
    % get distractor starts
    f=alignment.movie_distractor>distractThresh;
    fstarts=find(diff(f)==1);
    fstarts=fstarts+1;
    temp=zeros(size(alignment.movie_distractor));
    temp(fstarts)=1;
    figure();
    plot(alignment.movie_distractor);
    hold all;
    plot(temp);
    legend({'movie distractor raw','get beginnings'});
    title('getting beginning of movie distractor');
    temp(isnan(alignment.movie_distractor))=nan;
    alignment.movie_distractor=temp;
    
    % Get time delay
    timeIncs=diff(alignment.timesfromarduino(alignment.timesfromarduino~=0));
    mo=mode(timeIncs);
    timeIncs(timeIncs==mo)=nan;
    bettermode=mode(timeIncs); % in ms
    bettermode=bettermode/1000; % in seconds
    
    % get cue ITI
    cueITI=mode(diff(find(alignment.cueZone_onVoff>cueThresh)));
    disp('This is cue ITI in seconds');
    disp(cueITI*bettermode);
    
    % find movie distractors where time delay from most recent 
    % LED distractor OR real cue was at least this cue ITI
    % minus jitterInSec
    jitterInInd=floor(jitterInSec/bettermode);
    threshITI=cueITI-jitterInInd;
    cueAndDistract=alignment.movie_distractor>distractThresh | alignment.cueZone_onVoff>cueThresh;
    cueAndDistractInds=find(cueAndDistract);
    distractInds=find(alignment.movie_distractor>distractThresh);
    ITIforDistracts=nan(1,length(distractInds));
    for i=1:length(distractInds)
        temp=distractInds(i)-cueAndDistractInds(cueAndDistractInds~=distractInds(i));
        temp(temp<=0)=max(distractInds);
        ITIforDistracts(i)=min(temp);
    end
    alignment.backup_movie_distractor=alignment.movie_distractor;
    alignment.movie_distractor=zeros(size(alignment.movie_distractor));
    alignment.movie_distractor(distractInds(ITIforDistracts>=threshITI))=1;
    alignment.movie_distractor(isnan(alignment.backup_movie_distractor))=nan;
    %alignment.movie_distractor(distractInds(ITIforDistracts<threshITI))=0;
    figure();
    plot(alignment.backup_movie_distractor,'Color','k');
    hold all;
    plot(alignment.cueZone_onVoff,'Color','k');
    plot(alignment.movie_distractor,'Color','r');
    legend({'movie distractor raw','cue Zone','only long ITI'});
    title('take only long ITI distractors');
end

settings=plotCueTriggered_settings();
settings.isOrchestra=1;
settings.plotfields={'movie_distractor', ...
                     'cueZone_onVoff',...
                     'pelletLoaded',...
                     'pelletPresented',...
                     'optoOn',...
                     'interlockSamples',...
                     'reachStarts',...
                     'reach_ongoing',...
                     'success_reachStarts',...
                     'drop_reachStarts',...
                     'miss_reachStarts',...
                     'pelletmissingreach_reachStarts',...
                     'eating'};
settings.plotevents={'optoZone',...
                     'isGrooming',...
                     'success_reachStarts_pawOnWheel',...
                     'drop_reachStarts_pawOnWheel',...
                     'miss_reachStarts_pawOnWheel',...
                     'pelletmissingreach_reachStarts',...
                     'pelletPresented',...
                     'success_reachStarts',...
                     'drop_reachStarts',...
                     'miss_reachStarts',...
                     'cueZone_onVoff',...
                     'movie_distractor'};
settings.histoplotfields={'movie_distractor',...
                          'pelletPresented',...
                          'reachStarts'};
settings.trialType_plotfields={'movie_distractor',...
                          'pelletPresented',...
                          'success_reachStarts',...
                          'drop_reachStarts',...
                          'miss_reachStarts',...
                          'reachStarts',...
                          'goodReaches'};

tbt=plotCueTriggeredBehavior(alignment,'movie_distractor',0,settings);
if onlyLongITI==1
    save([alignmentDir '/tbt_distractor_onlyLongITI.mat'],'tbt');
%     save([alignmentDir '/plotSettings_distractor.mat'],'settings');
else
    save([alignmentDir '/tbt_distractor.mat'],'tbt');
    save([alignmentDir '/plotSettings_distractor.mat'],'settings');
end