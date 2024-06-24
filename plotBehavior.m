function tbt=plotBehavior(varargin)

if length(varargin)>4
    tbt=varargin{1};
    nameOfCue=varargin{2};
    excludePawOnWheelTrials=varargin{3};
    useTheseTrials=varargin{4};
    settings=varargin{5};
else
    tbt=varargin{1};
    nameOfCue=varargin{2};
    excludePawOnWheelTrials=varargin{3};
    useTheseTrials=varargin{4};
    % Get settings for this analysis
    settings=plotCueTriggered_settings();
end

% nameOfCue should be 'cue' for real cue
% 'cueZone_onVoff' for cue from movie
% or 'arduino_distractor' for distractor

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

plot_cues=plot_cues(ismember(plot_cues,useTheseTrials));

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
if ismember('pelletLoaded',fieldnames(tbt))
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
    endoftrialind=find(all(~isnan(tbt.(nameOfCue)),1)==1,1,'last');
end

% Also plot experiment as events in a scatter plot
figure();
k=1;
plotfields=settings.plotevents;
lastTrialShaded=0;
trialTypes=nan(1,length(plot_cues));
for i=plot_cues
    % Classify this trial type
    if (k==1 && settings.excludeFirstTrial==1) || (lastTrialShaded==0)
        % Skip this trial
    else
        % Classify trial type based on opto and shading of last trial
        if any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==1
            % opto and last trial was slow block, i.e., licking
            trialTypes(k)=1; % licking opto
        elseif ~any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==1
            % control and last trial was slow block, i.e., licking
            trialTypes(k)=2; % licking control
        elseif any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==2
            % opto and last trial was fast block, i.e., reaching
            trialTypes(k)=3;
        elseif ~any(tbt.optoOn(i,1:endoftrialind)>0.5) && lastTrialShaded==2
            % control and last trial was fast block, i.e., reaching
            trialTypes(k)=4;
        end
    end
    if ~isempty(settings.shading_type)
        % Shade some trials
        if ismember('ITI',settings.shading_type)
            % Shade trials according to ITI lengths
            shades=settings.shading_colors{find(ismember(settings.shading_type,'ITI'))};
            event_thresh=0.5;
            temp=tbt.(nameOfCue);
            event_ind_cue=find(temp(i,:)>event_thresh,1,'first');
            if isfield(tbt,'pelletPresented')
                event_ind_pellet=find(tbt.pelletPresented(i,:)>event_thresh);
            else
                event_ind_pellet=event_ind_cue;
            end
            if isempty(event_ind_pellet) || isempty(event_ind_cue)
            elseif event_ind_pellet(end)>length(timespertrial) || event_ind_cue>length(timespertrial)
            elseif any((timespertrial(event_ind_pellet)-timespertrial(event_ind_cue))>0 & (timespertrial(event_ind_pellet)-timespertrial(event_ind_cue))<settings.blockITIThresh)
                % Fast block
                lastTrialShaded=2;
                if strcmp(shades{2},'none')
                else
                    line([0 timespertrial(end)],[k k],'Color',shades{2},'LineWidth',10);
                end
            else
                % Slow block
                lastTrialShaded=1;
                if strcmp(shades{1},'none')
                else
                    line([0 timespertrial(end)],[k k],'Color',shades{1},'LineWidth',10);
                end
            end
%             elseif (timespertrial(event_ind_pellet(end))-timespertrial(event_ind_cue))>settings.blockITIThresh
%                 % Slow block
%                 if strcmp(shades{1},'none')
%                 else
%                     line([0 timespertrial(end)],[k k],'Color',shades{1},'LineWidth',10);
%                 end
%             else
%                 % Fast block
%                 if strcmp(shades{2},'none')
%                 else
%                     line([0 timespertrial(end)],[k k],'Color',shades{2},'LineWidth',10);
%                 end
%             end
        end
    end
    for j=1:length(plotfields)
        if ~isfield(tbt,plotfields{j})
            error([plotfields{j} ' field absent from tbt. See plotCueTriggered_settings.m to specify fields to plot.']);
        end
        currEvents=tbt.(plotfields{j});
        event_thresh=settings.eventThresh{j};
        event_ind=find(currEvents(i,:)>event_thresh);
        n=length(event_ind);
        if ischar(settings.firstN{j})
            if strcmp('all',settings.firstN{j})
                % plot all events
                if n>500
                    event_ind=event_ind(1:10:end);
                    n=length(event_ind);
                end
            end
        else
            % plot first n events
            n=settings.firstN{j};
        end
        if isempty(event_ind)
            continue
        end
        if n>length(event_ind)
            n=length(event_ind);
        end
        for l=1:n
            scatter([timespertrial(event_ind(l))-settings.shiftBack{j} timespertrial(event_ind(l))-settings.shiftBack{j}],[k k],[],'MarkerEdgeColor',settings.eventOutlines{j},...
                'MarkerFaceColor',settings.eventColors{j},...
                'LineWidth',settings.eventLineWidth);
            hold on;
        end       
    end
    k=k+1;
end

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

if settings.useFixedTimeBins==1
    % tbt.optoOn=double(tbt.optoZone>settings.eventThresh{1});
    disp('Started resampling for fixed time bins');
    % Put tbt into fixed time bins
    f=fieldnames(tbt);
    binTimes=settings.binMids;
    for i=1:length(f)
        disp(['Resampling ' f{i}]);
        if islogical(tbt.(f{i}))
            tbt.(f{i})=single(tbt.(f{i}));
        end
        temp=tbt.(f{i});
        newTimeTbt_temp=nan(size(temp,1),length(binTimes));
        for j=1:size(temp,1)
            temp2=temp(j,:);
            if isnan(temp2(1))
                temp2(~isnan(tbt.times(j,:)))=0;
            end
            % cut off anything that does not match times
            temp2(isnan(tbt.times(j,:)))=nan;
            curr=timeseries(temp2(~isnan(temp2)),tbt.times(j,~isnan(temp2)));
            if length(curr)<2
                res_curr.data=nan(1,size(newTimeTbt_temp,2));
            else
                res_curr=resample(curr,binTimes);
            end
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

function settings=plotCueTriggered_settings()

% Inter-trial interval (ITI) settings for this experiment
% Will help with choosing one cue timepoint per trial
% in case of issues with aliasing of instantaneous cue
settings.maxITI=50; % in seconds, maximal ITI
settings.minITI=2; % in seconds, minimal ITI

% Min height of cue or pellet presented peak relative to neighboring
% samples
% Increase for more stringent cue / pellet presented detection
% Decrease to pick up more peaks
% settings.relativePeakHeight=1*10^35; 
settings.nStdDevs=1; % relative peak height in terms of standard deviations away from mean
% i.e., greater than nStdDevs * nanstd(data) away from mean

% Experiment block types:

% Get to choose one or more types of trial shading for the event plot
% shading_type can include
% 'ITI'     blocks determined by different ITI lengths
settings.shading_type={'ITI'};
% settings.shading_type={};

% For experiments with blocks of different trial lengths
% If 'ITI' is included in shading_type, will shade each block of a different
% trial length with a different  color
settings.blockITIThresh=[0]; % boundaries between ITI block lengths, in seconds
% e.g., if there are 2 trial lengths, blockITIThresh is the threshold that
% distinguishes short and long trials

% Shading colors
% The indices here correspond to indices into shading_type
settings.shading_colors={{[0.9 0.9 0.9],'none'}};                        

% How many time points to include from previous trial at the beginning of
% each trial
settings.pointsFromPreviousTrial=100;


% Plotting trial-by-trial averages:

% Which fields from data to plot as averaged trial-by-trial
% List the field names to plot (these are field names from data variable in
% plotCueTriggeredBehavior.m)
settings.plotfields={'cueZone_onVoff', ...
                     'arduino_distractor',...
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
% settings.plotfields={'cueZone_onVoff', ...
%                      'arduino_distractor',...
%                      'pelletLoaded',...
%                      'pelletPresented',...
%                      'reachStarts',...
%                      'reach_ongoing',...
%                      'success_reachStarts',...
%                      'drop_reachStarts',...
%                      'miss_reachStarts',...
%                      'pelletmissingreach_reachStarts',...
%                      'eating'};

                 
% Plotting events:

% Which fields from data to plot as events
% Include field names in this list
% Will plot these events in this order, so event types later in the list
% will be plotted on top of the event types earlier in the list
settings.plotevents={'optoOn',...
                     'isGrooming',...
                     'success_reachStarts_pawOnWheel',...
                     'drop_reachStarts_pawOnWheel',...
                     'miss_reachStarts_pawOnWheel',...
                     'pelletmissingreach_reachStarts',...
                     'pelletPresented',...
                     'success_reachStarts',...
                     'drop_reachStarts',...
                     'miss_reachStarts',...
                     'isHold',...
                     'pawOnWheel',...
                     'cueZone_onVoff'};
% settings.plotevents={'success_reachStarts_pawOnWheel',...
%                      'drop_reachStarts_pawOnWheel',...
%                      'miss_reachStarts_pawOnWheel',...
%                      'pelletmissingreach_reachStarts',...
%                      'pelletPresented',...
%                      'success_reachStarts',...
%                      'drop_reachStarts',...
%                      'miss_reachStarts',...
%                      'cueZone_onVoff'};
                 
% Color for each event type
pawwheel_color='y'; % will fill in circle with this color if paw started from wheel for this reach
% Indices into eventColors and eventOutlines correspond to indices into plotevents
% eventColors gives fill of circle for each event type
% eventOutlines gives outline of circle for each event type
settings.eventColors={  'none',...
                        'none',...
                        pawwheel_color,...
                        pawwheel_color,...
                        pawwheel_color...
                        [0.8 0.8 0.8],...
                        'k',...
                        [0 0.75 0],...
                        'r',...
                        'c',...
                        'none',...
                        [0.8 0.8 0.8],...
                        'b'};
settings.eventOutlines={'m',...
                        [0.65 0.65 0.65],...
                        [0 0.75 0],...
                        'r',...
                        'c',...
                        [0.8 0.8 0.8],...
                        'none',...
                        [0 0.75 0],...
                        'none',...
                        'none',...
                        pawwheel_color,...
                        pawwheel_color,...
                        'b'};
% settings.eventColors={  pawwheel_color,...
%                         pawwheel_color,...
%                         pawwheel_color...
%                         [0.8 0.8 0.8],...
%                         'k',...
%                         [0 0.75 0],...
%                         'r',...
%                         'c',...
%                         'b'};
% settings.eventOutlines={[0 0.75 0],...
%                         'r',...
%                         'c',...
%                         [0.8 0.8 0.8],...
%                         'none',...
%                         [0 0.75 0],...
%                         'none',...
%                         'none',...
%                         'b'};
                    
% Event thresholds for each event type
% All indices with data values above the threshold will be counted as events
% May want to specify different thresholds for each event type
% Indices into eventThresh correspond to indices into plotevents
settings.eventThresh={0.5,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.2,...
                      0.2,...
                      0.2,...
                      0.5,...
                      0.5,...
                      0.5};
% settings.eventThresh={0.5,...
%                       0.5,...
%                       0.5,...
%                       0.5,...
%                       0.5,...
%                       0.2,...
%                       0.2,...
%                       0.2,...
%                       0.5};

% Will only plot firstN events of each type in each trial
% For example, will only plot first moment of cue on in each trial
% Indices into firstN correspond to indices into plotevents
% If value is 'all', will plot all events of this type in trial
settings.firstN={'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 1};
% settings.firstN={'all',...
%                  'all',...
%                  'all',...
%                  'all',...
%                  'all',...
%                  'all',...
%                  'all',...
%                  'all',...
%                  1};

% Width of outline of event marker
settings.eventLineWidth=1.5;

% Some events, like the cue, are detected as the peak of the resampled cue
% data. The onset of the cue may precede this peak by a fixed time. All
% events will be shifted backwards in time by shiftBack seconds, in the figure. 
% For example, if shiftBack is 0.2 for plotevents 'cue', the cue event will be
% plotted 200 ms earlier. Indices into shiftBack correspond to indices 
% into plotevents. If in doubt, set shiftBack to zero.
settings.shiftBack={0,...
                    0,...
                    0,...
                    0,...
                    0,...
                    0,...
                    0.2,...
                    0,...
                    0,...
                    0,...
                    0,...
                    0,...
                    0};
% settings.shiftBack={0,...
%                     0,...
%                     0,...
%                     0,...
%                     0.2,...
%                     0,...
%                     0,...
%                     0,...
%                     0.2};
                
% Plot histogram of events (or trial-by-trial average), with data types
% overlapping:
                
% Which fields from data to plot as overlapping averaged trial-by-trial 
% List the field names to plot (these are field names from data variable in
% plotCueTriggeredBehavior.m) in histoplotfields
% Note that will also shift each event time for each data type by 
% histoshiftBack
settings.histoDropGrooming=0; % if 1, will drop all times periods in which mouse was grooming from histograms, otherwise 0
settings.dropChewingInCue=0; % if 1, will drop all trials when mouse was already chewing during cue
settings.histoplotfields={'cueZone_onVoff',...
                          'pelletPresented',...
                          'reachStarts'};
settings.histoshiftBack={0,...
                         0.2,...
                         0};


% For comparing successes and failures with different trial types
settings.excludeFirstTrial=0; % if 1, will drop first trial, which is usually incomplete
settings.binByN=6; % how much to down-sample tbt in time
settings.trialType_plotfields={'cueZone_onVoff',...
                          'pelletPresented',...
                          'success_reachStarts',...
                          'drop_reachStarts',...
                          'miss_reachStarts',...
                          'reachStarts',...
                          'goodReaches'};
settings.trialType_shiftBack={0,...
                         0.2,...
                         0,...
                         0,...
                         0,...
                         0,...
                         0};
settings.trialType_name={'licking opto',...
                         'licking control',...
                         'opto reaching',...
                         'control reaching'};
                     
% Whether to use fixed time bins for histograms, or just use time bins emerging
% from previous analysis
settings.useFixedTimeBins=0; % 1 if yes, 0 if no
settings.binMids=0:0.035:17; % in seconds


end
