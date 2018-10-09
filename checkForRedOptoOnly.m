function wasControl=checkForRedOptoOnly(out)

event_thresh=0.5;
perc_thresh=10; % would not have set cue missing probability to less than this %

trials_cueOn=any(out.cueOn>event_thresh,2);
trials_cueOn=single(trials_cueOn);

trials_optoOn=any(out.optoOn>event_thresh,2);
trials_optoOn=single(trials_optoOn);

optoWithoutCue=trials_optoOn==1 & trials_cueOn==0;
disp('Fraction of trials with opto on but cue did not go on');
frac=sum(optoWithoutCue)./length(optoWithoutCue);
disp(frac);
if frac*100<perc_thresh
    wasControl=0;
    disp('This session WAS NOT opto ONLY control');
else
    wasControl=1;
    disp('This session WAS opto ONLY control');
end