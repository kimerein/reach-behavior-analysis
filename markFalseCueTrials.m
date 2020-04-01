function alignment=markFalseCueTrials(alignment)

% enters false cues into field 'falseCueOn'

% find false cues
% note that, normally, pellets are loaded and then presented to mouse 5
% trials later (cue should follow load by 5 trials)
% if cue is absent on this 5th trial, then doing "no cue" control

event_thresh=0.5;
alignment.falseCueOn=zeros(size(alignment.cue));
pelletPresented=alignment.pelletPresented;
pelletLoaded=alignment.pelletLoaded;
cue=alignment.cue;
[~,loads]=findpeaks(pelletLoaded,'MinPeakProminence',event_thresh);
[~,presents]=findpeaks(pelletPresented,'MinPeakProminence',event_thresh);
[~,cues]=findpeaks(cue,'MinPeakProminence',event_thresh);
% find presents just preceding cues
presents_preceding=nan(1,length(cues));
for i=1:length(cues)
    currcue=cues(i);
    temp=currcue-presents;
    temp(temp<0)=nan;
    [~,mind]=nanmin(temp);
    presents_preceding(i)=presents(mind);
end
presentToCueDelay=mode(cues-presents_preceding);
% for trials where pellet was presented but cue did not go off, add false
% cue
didload=nan(1,length(presents));
for i=1:length(presents)-5
    % iterate through each wheel turn
    currpresent=presents(i);
    nextpresent=presents(i+1);
    % did load this trial?
    if any(loads>currpresent & loads<nextpresent)
        didload(i)=1;
    else
        didload(i)=0;
    end
end
didcue=nan(1,length(presents));
for i=6:length(presents)
    % iterate through each wheel turn
    currpresent=presents(i);
    if i+1>length(presents)
        nextpresent=length(pelletPresented);
    else
        nextpresent=presents(i+1);
    end
    if any(cues>currpresent & cues<nextpresent)
        didcue(i)=1;
    else
        didcue(i)=0;
    end
end
trialsWithPelletWithoutCue=didcue(6:end)-didload(1:end-5);
% mark these trials as false cues
temp=6:length(presents);
presentIndsForFalseCues=temp(trialsWithPelletWithoutCue~=0);
presentsForFalseCues=presents(presentIndsForFalseCues);
alignment.falseCueOn(presentsForFalseCues+presentToCueDelay)=1;
figure();
plot(alignment.pelletPresented,'Color','k'); 
hold on; 
plot(alignment.cue,'Color','b');
plot(alignment.falseCueOn,'Color','r');
leg={'pellet presented','cue','false cue'};
xlabel('indices');
ylabel('events');
legend(leg);
title('detecting false cues');


