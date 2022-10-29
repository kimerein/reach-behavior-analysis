function aligned=manuallyAddInCuesGUI(aligned)

figure();
plot(aligned.cueZone,'Color','k');
hold on;
plot(aligned.cue*range(aligned.cueZone)+min(aligned.cueZone),'Color','r');
plot(aligned.cueZone_onVoff*range(aligned.cueZone)+min(aligned.cueZone),'Color','b');
legend({'original from movie','from arduino','fixed'});

addCues=input('Add indices of cues to add as array, e.g., [x y z]: ');
aligned.cueZone_onVoff(addCues)=nanmax(aligned.cueZone_onVoff);

figure();
plot(aligned.cueZone,'Color','k');
hold on;
plot(aligned.cue*range(aligned.cueZone)+min(aligned.cueZone),'Color','r');
plot(aligned.cueZone_onVoff*range(aligned.cueZone)+min(aligned.cueZone),'Color','b');
legend({'original from movie','from arduino','with cues added'});

end

