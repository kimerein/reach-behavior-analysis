function tbt=plotCueTriggeredBehavior(varargin)

if length(varargin)>3
    data=varargin{1};
    nameOfCue=varargin{2};
    excludePawOnWheelTrials=varargin{3};
    settings=varargin{4};
else
    data=varargin{1};
    nameOfCue=varargin{2};
    excludePawOnWheelTrials=varargin{3};
    % Get settings for this analysis
    settings=plotCueTriggered_settings();
end

% nameOfCue should be 'cue' for real cue
% 'cueZone_onVoff' for cue from movie
% or 'arduino_distractor' for distractor

% Make movie distractor only the START of movie distractor -- change to
% instantaneous
if strcmp(nameOfCue,'movie_distractor')
    f=data.movie_distractor>0.5;
    fstarts=find(diff(f)==1);
    fstarts=fstarts+1;
    temp=zeros(size(data.movie_distractor));
    temp(fstarts)=1;
    figure();
    plot(data.movie_distractor);
    hold all;
    plot(temp);
    legend({'movie distractor raw','get beginnings'});
    title('getting beginning of movie distractor');
    temp(isnan(data.movie_distractor))=nan;
    data.movie_distractor=temp;
elseif strcmp(nameOfCue,'optoOnly')
    f=data.optoOnly>0.5;
    fstarts=find(diff(f)==1);
    fstarts=fstarts+1;
    temp=zeros(size(data.optoOnly));
    temp(fstarts)=1;
    figure();
    plot(data.optoOnly);
    hold all;
    plot(temp);
    legend({'optoOnly raw','get beginnings'});
    title('getting beginning of optoOnly');
    temp(isnan(data.optoOnly))=nan;
    data.optoOnly=temp;
end

% If running automatically on Harvard server, do not use manual optoZone
% threshold
if settings.isOrchestra==1 & isfield(data,'optoZone')
    threshForOnVsOff=min(data.optoZone)+settings.fractionRange*range(data.optoZone);
    for i=1:length(settings.plotevents)
        if strcmp(settings.plotevents{i},'optoZone')
            settings.eventThresh{i}=threshForOnVsOff;
            break
        end
    end
end

% Get cue/data type for triggering trial-by-trial data
cue=data.(nameOfCue);

% In case of issues with aliasing of instantaneous cue
maxITI=settings.maxITI; % in seconds, maximal ITI
minITI=settings.minITI; % in seconds, minimal ITI

% Get time delay
timeIncs=diff(data.timesfromarduino(data.timesfromarduino~=0));
mo=mode(timeIncs);
timeIncs(timeIncs==mo)=nan;
bettermode=mode(timeIncs); % in ms
% bettermode=mo;
bettermode=bettermode/1000; % in seconds

% Fix aliasing issues with resampled data
% if strcmp(nameOfCue,'cueZone_onVoff')
if strcmp(nameOfCue,'cue') | strcmp(nameOfCue,'cueZone_onVoff') | strcmp(nameOfCue,'falseCueOn') | strcmp(nameOfCue,'movie_distractor') | strcmp(nameOfCue,'optoOnly')
    [cue,cueInds,cueIndITIs]=fixAlias_forThreshCue(cue,maxITI,minITI,bettermode);
else
    [cue,cueInds,cueIndITIs]=fixAliasing(cue,maxITI,minITI,bettermode);
end
[data.pelletPresented,presentedInds]=fixAliasing(data.pelletPresented,maxITI,minITI,bettermode);

% Checking cue detection
figure();
plot(cue./nanmax(cue));
hold on;
plot(data.pelletPresented./nanmax(data.pelletPresented),'Color','k');
for i=1:length(cueInds)
    scatter(cueInds(i),1,[],'r');
end
for i=1:length(presentedInds)
    scatter(presentedInds(i),1,[],[0.5 0.5 0.5]);
end
title('Checking cue selection');
legend({'Cue','Pellet presented','Cue detected','Pellet presented detected'});

data.timesfromarduino=data.timesfromarduino./1000; % convert from ms to seconds

% Turn pellet presented into events
pelletPresented=data.pelletPresented;
sigmax=max(max(pelletPresented));
sigthresh=sigmax/2;
temp=zeros(size(pelletPresented));
temp(pelletPresented>=sigthresh)=1;
data.pelletPresented=temp;

% Set up trial-by-trial, tbt
pointsFromPreviousTrial=settings.pointsFromPreviousTrial;
f=fieldnames(data);
for i=1:length(f)
    tbt.(f{i})=nan(length(cueInds),max(cueIndITIs)+pointsFromPreviousTrial);
end
% Add times
tbt.times=nan(length(cueInds),max(cueIndITIs)+pointsFromPreviousTrial);

% Organize trial-by-trial data
for i=1:length(cueInds)
    if i==length(cueInds)
        theseInds=cueInds(i)-pointsFromPreviousTrial:length(cue);
    elseif i==1
        theseInds=1:cueInds(i+1)-1;
    else
        if cueInds(i)-pointsFromPreviousTrial<0
            theseInds=1:cueInds(i+1)-1;
        else
            theseInds=cueInds(i)-pointsFromPreviousTrial:cueInds(i+1)-1;
        end
    end
    for j=1:length(f)
        temp=tbt.(f{j});
        temp1=data.(f{j});
        temp(i,1:length(theseInds))=temp1(theseInds);
        tbt.(f{j})=temp;
    end
    % Get times
    tbt.times(i,:)=nan;
    tbt.times(i,1:length(theseInds))=data.movieframeinds(theseInds).*bettermode;
end

% Zero out nans
% for i=1:length(f)
%     temp=tbt.(f{i});
%     temp(isnan(temp))=0;
%     tbt.(f{i})=temp;
% end

% Get times per trial
tbt.times=tbt.times-repmat(nanmin(tbt.times,[],2),1,size(tbt.times,2));
timespertrial=nanmean(tbt.times,1);

% Exclude trials where paw was on wheel while wheel turning
if excludePawOnWheelTrials==1
    % Find trials where paw was on wheel while wheel turning
    plot_cues=[];
    for i=1:size(tbt.(nameOfCue),1)
        presentInd=find(tbt.pelletPresented(i,:)>0.5,1,'first');
        temp=tbt.(nameOfCue);
        cueInd=find(temp(i,:)>0.5,1,'first');
        pawWasOnWheel=0;
%         if any(tbt.pawOnWheel(i,presentInd:cueInd)>0.5)
        if any(tbt.reach_ongoing(i,presentInd:cueInd)>0.5) || any(tbt.isHold(i,presentInd:cueInd)>0.5)
            pawWasOnWheel=1;
        else
            plot_cues=[plot_cues i];
        end
    end
else
    plot_cues=1:size(tbt.(nameOfCue),1);
end
if settings.excludeFirstTrial==1
    plot_cues=plot_cues(~ismember(plot_cues,1));
end

% Plot trial-by-trial average
figure();
plotfields=settings.plotfields;
ha=tight_subplot(length(plotfields),1,[0.06 0.03],[0.05 0.05],[0.1 0.03]);
for i=1:length(plotfields)
    currha=ha(i);
    axes(currha);
    if ~isfield(tbt,plotfields{i})
        error([plotfields{i} ' field absent from tbt. See plotCueTriggered_settings.m to specify fields to plot.']);
    end
    temp=tbt.(plotfields{i});
    plot(timespertrial,nanmean(temp(plot_cues,:),1));
    title(plotfields{i},'Interpreter','none');
end

% Ask user when min trial ends based on distribution of wheel turn times
if settings.isOrchestra~=1
    figure(); plot(nanmean(tbt.pelletLoaded,1));
    hold on; plot(nanmean(tbt.cueZone_onVoff,1),'Color','b');
    hold on; plot(nanmean(tbt.pelletPresented,1),'Color','k');
    leg={'pellet loaded','cue','pellet presented'};
    xlabel('indices');
    ylabel('av');
    title('average across trials');
    legend(leg);
    endoftrialind=input('Enter index showing minimum trial length -- note that any opto within a trial should occur before this index. ');
    if ~isnumeric(endoftrialind)
        error('Please enter a number.');
    end
else
    endoftrialind=250;
end

% Also plot experiment as events in a scatter plot
figure();
k=1;
plotfields=settings.plotevents;
lastTrialShaded=0;
trialTypes=nan(1,length(plot_cues));
% for i=plot_cues
%     % Classify this trial type
%     if (k==1 && settings.excludeFirstTrial==1) || (lastTrialShaded==0)
%         % Skip this trial
%     else
%         % Classify trial type based on opto and shading of last trial
%         if any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==1
%             % opto and last trial was slow block, i.e., licking
%             trialTypes(k)=1; % licking opto
%         elseif ~any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==1
%             % control and last trial was slow block, i.e., licking
%             trialTypes(k)=2; % licking control
%         elseif any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==2
%             % opto and last trial was fast block, i.e., reaching
%             trialTypes(k)=3;
%         elseif ~any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==2
%             % control and last trial was fast block, i.e., reaching
%             trialTypes(k)=4;
%         end
%     end
%     if ~isempty(settings.shading_type)
%         % Shade some trials
%         if ismember('ITI',settings.shading_type)
%             % Shade trials according to ITI lengths
%             shades=settings.shading_colors{find(ismember(settings.shading_type,'ITI'))};
%             event_thresh=0.5;
%             temp=tbt.(nameOfCue);
%             event_ind_cue=find(temp(i,:)>event_thresh,1,'first');
%             event_ind_pellet=find(tbt.pelletPresented(i,:)>event_thresh);
%             if isempty(event_ind_pellet) || isempty(event_ind_cue)
%             elseif event_ind_pellet(end)>length(timespertrial) || event_ind_cue>length(timespertrial)
%             elseif any((timespertrial(event_ind_pellet)-timespertrial(event_ind_cue))>0 & (timespertrial(event_ind_pellet)-timespertrial(event_ind_cue))<settings.blockITIThresh)
%                 % Fast block
%                 lastTrialShaded=2;
%                 if strcmp(shades{2},'none')
%                 else
%                     line([0 timespertrial(end)],[k k],'Color',shades{2},'LineWidth',10);
%                 end
%             else
%                 % Slow block
%                 lastTrialShaded=1;
%                 if strcmp(shades{1},'none')
%                 else
%                     line([0 timespertrial(end)],[k k],'Color',shades{1},'LineWidth',10);
%                 end
%             end
% %             elseif (timespertrial(event_ind_pellet(end))-timespertrial(event_ind_cue))>settings.blockITIThresh
% %                 % Slow block
% %                 if strcmp(shades{1},'none')
% %                 else
% %                     line([0 timespertrial(end)],[k k],'Color',shades{1},'LineWidth',10);
% %                 end
% %             else
% %                 % Fast block
% %                 if strcmp(shades{2},'none')
% %                 else
% %                     line([0 timespertrial(end)],[k k],'Color',shades{2},'LineWidth',10);
% %                 end
% %             end
%         end
%     end
%     for j=1:length(plotfields)
%         if ~isfield(tbt,plotfields{j})
%             error([plotfields{j} ' field absent from tbt. See plotCueTriggered_settings.m to specify fields to plot.']);
%         end
%         currEvents=tbt.(plotfields{j});
%         event_thresh=settings.eventThresh{j};
%         event_ind=find(currEvents(i,:)>event_thresh);
%         n=length(event_ind);
%         if ischar(settings.firstN{j})
%             if strcmp('all',settings.firstN{j})
%                 % plot all events
%                 if n>500
%                     event_ind=event_ind(1:10:end);
%                     n=length(event_ind);
%                 end
%             end
%         else
%             % plot first n events
%             n=settings.firstN{j};
%         end
%         if isempty(event_ind)
%             continue
%         end
%         for l=1:n
%             scatter([timespertrial(event_ind(l))-settings.shiftBack{j} timespertrial(event_ind(l))-settings.shiftBack{j}],[k k],[],'MarkerEdgeColor',settings.eventOutlines{j},...
%                 'MarkerFaceColor',settings.eventColors{j},...
%                 'LineWidth',settings.eventLineWidth);
%             hold on;
%         end       
%     end
%     k=k+1;
% end

% Drop grooming time periods?
groomingTrials=[];
% Drop trials when mouse was chewing during cue?
chewingTrials=[];
if settings.histoDropGrooming==1 || settings.dropChewingInCue==1
    % Find trials when animal was grooming during cue and exclude these
    % trials
    % Find trials when mouse was already chewing when cue came on
    % Exclude these trials
    temp=tbt.(nameOfCue);
    for i=1:size(temp,1)
        cueInd=find(temp(i,:)>0.5,1,'first');
        if tbt.isGrooming(i,cueInd)>0.5 % mouse is grooming during cue
            % exclude this trial
            groomingTrials=[groomingTrials i];
        end
        if tbt.isChewing(i,cueInd)>0.5 % mouse is chewing during cue
%         if tbt.isChewing(i,cueInd-5)>0.5 % mouse is chewing during cue
            % exclude this trial
            chewingTrials=[chewingTrials i];
        end
    end
    if settings.histoDropGrooming==1
        plot_cues=plot_cues(~ismember(plot_cues,groomingTrials));
        disp('Excluding from histograms the following trials where mouse was grooming during cue');
        disp(groomingTrials);
    end
    if settings.dropChewingInCue==1
        plot_cues=plot_cues(~ismember(plot_cues,chewingTrials));
        disp('Excluding from histograms the following trials where mouse was chewing during cue');
        disp(chewingTrials);
    end
end

% Plot overlap trial-by-trial average
figure();
plotfields=settings.histoplotfields;
for i=1:length(plotfields)
    temp=tbt.(plotfields{i});
    if i==1
        plot(timespertrial-settings.histoshiftBack{i},nansum(temp(plot_cues,:),1));
    else
        temp2=nansum(temp(plot_cues,:),1);
        plot(timespertrial-settings.histoshiftBack{i},temp2.*(ma/nanmax(temp2)));
    end
    hold all;
    if i==1
        ma=nanmax(nansum(temp(plot_cues,:),1));
    end
end
legend(plotfields);

% Fill in times with respect to cue in each trial (with respect to trial
% start)
tbt.times_wrt_trial_start=nan(size(tbt.times));
for i=1:size(tbt.times)
    useTimesForSeries=tbt.times(i,:);
    temp=tbt.(nameOfCue);
    cueInd=find(temp(i,:)>0.5,1,'first');
    useTimesForSeries=useTimesForSeries-useTimesForSeries(cueInd); % cue is now occurring at time zero
    useTimesForSeries=useTimesForSeries+pointsFromPreviousTrial*bettermode; % adjust zero to be start of trial
    % times after monotonically increasing times are no longer in
    % this trial
    isMonotonicIncrease=([useTimesForSeries(2:end) 0]-[useTimesForSeries(1:end-1) 0])>0;
    useTimesForSeries(~isMonotonicIncrease)=nan;
    tbt.times_wrt_trial_start(i,:)=useTimesForSeries;
end

if settings.useFixedTimeBins==1
    % tbt.optoOn=double(tbt.optoZone>settings.eventThresh{1});
    disp('Started resampling for fixed time bins');
    % Put tbt into fixed time bins
    f=fieldnames(tbt);
    binTimes=settings.binMids;
    for i=1:length(f)
        disp(['Resampling ' f{i}]);
        temp=tbt.(f{i});
        newTimeTbt_temp=nan(size(temp,1),length(binTimes));
        for j=1:size(temp,1)
            temp2=temp(j,:);
            if isnan(temp2(1))
                temp2(~isnan(tbt.times(j,:)))=0;
            end
            % cut off anything that does not match times
            temp2(isnan(tbt.times(j,:)))=nan;
            useTimesForSeries=tbt.times_wrt_trial_start(j,:);
            curr=timeseries(temp2(~isnan(temp2) & ~isnan(useTimesForSeries)),useTimesForSeries(~isnan(temp2) & ~isnan(useTimesForSeries)));
            res_curr=resample(curr,binTimes);
            newTimeTbt_temp(j,1:length(res_curr.data))=res_curr.data;
        end
        newTimeTbt.(f{i})=newTimeTbt_temp;
    end
    disp('Finished resampling for fixed time bins');
    tbt=newTimeTbt;
else
    % Time down-sample tbt
    f=fieldnames(tbt);
    binWins=1:settings.binByN:size(tbt.(f{1}),2);
    temp2=nan(size(tbt.(f{1}),1),length(binWins)-1);
    for i=1:length(f)
        temp=tbt.(f{i});
        for j=1:length(binWins)-1
            temp2(:,j)=sum(temp(:,binWins(j):binWins(j+1)-1),2);
        end
        ds_tbt.(f{i})=temp2;
    end
    backup_tbt=tbt;
    tbt=ds_tbt;
    new_timespertrial=nan(1,length(binWins)-1);
    for j=1:length(binWins)-1
        new_timespertrial(j)=nanmean(timespertrial(binWins(j):binWins(j+1)-1));
    end
    timespertrial=new_timespertrial;
    
    % Add together reaches when pellet available, starting from perch
    tbt.goodReaches=tbt.success_reachStarts+tbt.drop_reachStarts+tbt.miss_reachStarts;
    
    u=unique(trialTypes);
    u=u(~isnan(u));
    backup_plot_cues=plot_cues;
    for j=1:length(u)
        plot_cues=backup_plot_cues(trialTypes==u(j)); % Only use trials of this type
        temp=1:size(tbt.(plotfields{1}),1);
        plot_cues=ismember(temp,plot_cues);
        % Plot overlap trial-by-trial average
        figure();
        plotfields=settings.trialType_plotfields;
        for i=1:length(plotfields)
            temp=tbt.(plotfields{i});
            if i==1
                plot(timespertrial-settings.trialType_shiftBack{i},nansum(temp(plot_cues,:),1)./sum(plot_cues));
            else
                temp2=nansum(temp(plot_cues,:),1);
                %             plot(timespertrial-settings.trialType_shiftBack{i},temp2.*(ma/nanmax(temp2)));
                plot(timespertrial-settings.trialType_shiftBack{i},temp2./sum(plot_cues));
            end
            hold all;
            if i==1
                ma=nanmax(nansum(temp(plot_cues,:),1));
            end
        end
        legend(plotfields);
        title(settings.trialType_name{j});
    end
    
    tbt=backup_tbt;
end

end

function [cue,cueInds,cueIndITIs]=fixAlias_forThreshCue(cue,maxITI,minITI,bettermode)

settings=plotCueTriggered_settings();
peakHeight=0.5;

[pks,locs]=findpeaks(cue);
cueInds=locs(pks>peakHeight);
% cueInds=[1 cueInds length(cue)]; % in case aliasing problem is at edges
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode>(maxITI*1.5));
for i=1:length(checkTheseIntervals)
    indsIntoCue=cueInds(checkTheseIntervals(i))+floor((maxITI/2)./bettermode):cueInds(checkTheseIntervals(i)+1)-floor((maxITI/2)./bettermode);
    if any(cue(indsIntoCue)>0.001)
        [~,ma]=max(cue(indsIntoCue));
        cue(indsIntoCue(ma))=max(cue);
    end
end

% [pks,locs]=findpeaks(cue);
[pks,locs]=findpeaks(cue,'MinPeakDistance',floor((minITI*0.75)/bettermode),'MinPeakProminence',peakHeight);
cueInds=locs(pks>peakHeight);
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode<(minITI*0.75));
if ~isempty(checkTheseIntervals)
    for i=1:length(checkTheseIntervals)
        cue(cueInds(checkTheseIntervals(i)))=0;
        cueInds(checkTheseIntervals(i))=nan;
    end
end
cueInds=cueInds(~isnan(cueInds));
cueIndITIs=diff(cueInds);

cue=cue./nanmax(cue);

end

function [cue,cueInds,cueIndITIs]=fixAliasing(cue,maxITI,minITI,bettermode)

cue=nonparamZscore(cue); % non-parametric Z score

settings=plotCueTriggered_settings();
peakHeight=nanmean(cue)+settings.nStdDevs*nanstd(cue);
relativePeakHeight=settings.nStdDevs*nanstd(cue);

[pks,locs]=findpeaks(cue);
cueInds=locs(pks>peakHeight);
% cueInds=[1 cueInds length(cue)]; % in case aliasing problem is at edges
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode>(maxITI*1.5));
for i=1:length(checkTheseIntervals)
    indsIntoCue=cueInds(checkTheseIntervals(i))+floor((maxITI/2)./bettermode):cueInds(checkTheseIntervals(i)+1)-floor((maxITI/2)./bettermode);
    if any(cue(indsIntoCue)>0.001)
        [~,ma]=max(cue(indsIntoCue));
        cue(indsIntoCue(ma))=max(cue);
    end
end

% [pks,locs]=findpeaks(cue);
[pks,locs]=findpeaks(cue,'MinPeakDistance',floor((minITI*0.75)/bettermode),'MinPeakProminence',relativePeakHeight);
cueInds=locs(pks>peakHeight);
cueIndITIs=diff(cueInds);
checkTheseIntervals=find(cueIndITIs*bettermode<(minITI*0.75));
if ~isempty(checkTheseIntervals)
    for i=1:length(checkTheseIntervals)
        cue(cueInds(checkTheseIntervals(i)))=0;
        cueInds(checkTheseIntervals(i))=nan;
    end
end
cueInds=cueInds(~isnan(cueInds));
cueIndITIs=diff(cueInds);

cue=cue./nanmax(cue);

end
