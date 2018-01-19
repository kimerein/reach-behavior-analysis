function plotCueTriggeredExpt(tbt,nameOfCue,excludePawOnWheelTrials)

% nameOfCue should be 'cue' for real cue
% 'cueZone_onVoff' for cue from movie
% or 'arduino_distractor' for distractor

% Get settings for this analysis
settings=plotCueTriggered_settings();

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
        if any(tbt.pawOnWheel(i,presentInd:cueInd)>0.5)
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
    temp=tbt.(plotfields{i});
    plot(timespertrial,nanmean(temp(plot_cues,:),1));
    title(plotfields{i},'Interpreter','none');
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
        if any(tbt.optoOn(i,:)>0.5) && lastTrialShaded==1
            % opto and last trial was slow block, i.e., licking
            trialTypes(k)=1; % licking opto
        elseif ~any(tbt.optoOn(i,:)>0.5) && lastTrialShaded==1
            % control and last trial was slow block, i.e., licking
            trialTypes(k)=2; % licking control
        elseif any(tbt.optoOn(i,:)>0.5) && lastTrialShaded==2
            % opto and last trial was fast block, i.e., reaching
            trialTypes(k)=3;
        elseif ~any(tbt.optoOn(i,:)>0.5) && lastTrialShaded==2
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
            event_ind_pellet=find(tbt.pelletPresented(i,:)>event_thresh);
            if isempty(event_ind_pellet)
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
        currEvents=tbt.(plotfields{j});
        event_thresh=settings.eventThresh{j};
        event_ind=find(currEvents(i,:)>event_thresh);
        n=length(event_ind);
        if ischar(settings.firstN{j})
            if strcmp('all',settings.firstN{j})
                % plot all events
            end
        else
            % plot first n events
            n=settings.firstN{j};
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
for j=1:length(u)
    plot_cues=trialTypes==u(j); % Only use trials of this type
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

end
