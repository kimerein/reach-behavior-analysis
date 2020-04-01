function aligned=cleanUpCue(aligned)

% Settings

% This is the fraction of range above the min to use as cut-off threshold
% for cue on
<<<<<<< HEAD
thresh=0.5; 
=======
% thresh=0.57; 
thresh=0.33; 
>>>>>>> origin/master

% If this is 1, subtract off LED distractor, because there is contamination
% of cue zone from LED distractor
subtractDistract=1;

cueZone=aligned.cueZone;
if subtractDistract==1
    mi=min(cueZone);
    cueZone=cueZone-(0.01*range(cueZone)*aligned.movie_distractor);
    cueZone(cueZone<mi)=mi;
end

temp=nan(size(cueZone));
isOn=cueZone>(thresh*range(cueZone))+min(cueZone);
temp(isOn)=1;
temp(~isOn)=0;
aligned.cueZone_onVoff=temp;

figure();
plot(cueZone,'Color','k'); 
hold on;
line([0 length(cueZone)],[(thresh*range(cueZone))+min(cueZone) (thresh*range(cueZone))+min(cueZone)],'Color','r');
title('Cleaning up cue from movie');