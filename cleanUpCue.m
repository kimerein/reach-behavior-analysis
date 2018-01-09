function aligned=cleanUpCue(aligned)

% This is the fraction of range above the min to use as cut-off threshold
% for cue on
thresh=0.75; 

temp=nan(size(aligned.cueZone));
isOn=aligned.cueZone>(thresh*range(aligned.cueZone))+min(aligned.cueZone);
temp(isOn)=1;
temp(~isOn)=0;
aligned.cueZone_onVoff=temp;