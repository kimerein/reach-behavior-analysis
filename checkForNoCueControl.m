function [wasControl,out]=checkForNoCueControl(out)

% out is from out=parseSerialOut_wrapper([videoFile(1:endofDir(end)) 'OUTPUT.txt'],[videoFile(1:endofVfname(end)-1) '_parsedOutput.mat'])
% returns whether or not this session was a "no cue" control

% note that, normally, pellets are loaded and then presented to mouse 5
% trials later (cue should follow load by 5 trials)
% if cue is absent on this 5th trial, then doing "no cue" control

event_thresh=0.5;
perc_thresh=3; % would not have set cue missing probability to less than this %

trials_pelletLoaded=any(out.pelletLoaded>event_thresh,2);
trials_cueOn=any(out.cueOn>event_thresh,2);

trials_pelletLoaded=single(trials_pelletLoaded);
trials_cueOn=single(trials_cueOn);

trialsWithPelletWithoutCue=trials_cueOn(6:end)-trials_pelletLoaded(1:end-5);
disp('Fraction of trials where pellet was presented but cue did not go on');
frac=sum(trialsWithPelletWithoutCue~=0)/length(trialsWithPelletWithoutCue);
disp(frac);
if frac*100<perc_thresh
    wasControl=0;
    disp('This session WAS NOT no cue control');
else 
    wasControl=1;
    disp('This session WAS no cue control');
end
out.event_thresh=event_thresh;
out.perc_thresh=perc_thresh;
out.perc_missing_cue=frac*100;



