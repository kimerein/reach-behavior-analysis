function alignment=correctFalseCues_basedOnPellet(alignment)

thresh=0.5; % thresh for detecting events
maxDelay=1; % in seconds, max time between pellet arrives and false cue on
% this is important, because sometimes pellet does not load

% Get time delay
timeIncs=diff(alignment.timesfromarduino(alignment.timesfromarduino~=0));
mo=mode(timeIncs);
timeIncs(timeIncs==mo)=nan;
bettermode=mode(timeIncs); % in ms
% bettermode=mo;
timestep=bettermode/1000; % in seconds
maxDelayInd=floor(maxDelay/timestep);

pelletPresent=alignment.pelletPresent;
falseCue=alignment.falseCueOn;

% find pellet arrives in movie just before each false cue
f=find(falseCue>thresh);
pelletArrivesBeforeFalseCue=nan(1,length(f));
for i=1:length(f)
    curr=f(i);
    % find time when pellet arrives just preceding this
    fi=find(pelletPresent(curr:-1:2)>thresh & pelletPresent(curr-1:-1:1)<thresh,1,'first');
    if isempty(fi)
        continue
    end
    pellet_f=curr-(fi-1);
    pelletArrivesBeforeFalseCue(i)=pellet_f;
end

% eliminate any trials for which pellet arrives more than maxDelay seconds before false cue
pelletArrivesBeforeFalseCue(f-pelletArrivesBeforeFalseCue>maxDelayInd)=nan;
usualDelay=mode(f-pelletArrivesBeforeFalseCue);
pelletArrivesBeforeFalseCue=pelletArrivesBeforeFalseCue(~isnan(pelletArrivesBeforeFalseCue));

alignment.falseCueOn_backup=alignment.falseCueOn;
alignment.pelletArrivesBeforeFalseCue=zeros(size(alignment.pelletPresent));
alignment.pelletArrivesBeforeFalseCue(pelletArrivesBeforeFalseCue)=1;
% change falseCueOn so that aligned with pellet arrival, if possible
alignment.falseCueOn=zeros(size(alignment.pelletPresent));
alignment.falseCueOn(pelletArrivesBeforeFalseCue+usualDelay)=1;

% Plot sanity check
figure();
for i=1:length(pelletArrivesBeforeFalseCue)
    plot(alignment.pelletPresent(pelletArrivesBeforeFalseCue(i)-5:pelletArrivesBeforeFalseCue(i)+5));
    hold all;
end
title('checking alignment of pellet arrival');


