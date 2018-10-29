function [aligned,settings]=cleanUpCue_basedOnArduino(aligned)
 
isOrchestra=1;

if isOrchestra==1
    minProm=prctile(aligned.cueZone,95)-prctile(aligned.cueZone,40);
    minProm2=prctile(aligned.cueZone,60)-prctile(aligned.cueZone,20);
else
    minProm=10000;
    minProm2=5000;
end
settings.minProm=minProm;
settings.minProm2=minProm2;

[~,fmovie,widths,proms]=findpeaks(aligned.cueZone,'MinPeakProminence',minProm);
[~,farduino]=findpeaks(aligned.cue,'MinPeakHeight',0.5);
usePeak=zeros(size(aligned.cueZone));
usePeak(isnan(aligned.cueZone))=nan;
for i=1:length(farduino)
    currArduinoPeak=farduino(i);
    [~,mi]=min(abs(currArduinoPeak-fmovie));
    for j=1:50
        temp=aligned.cueZone(fmovie(mi)-(j-1))-aligned.cueZone(fmovie(mi)-j);
        if temp>minProm2
            break
        end
    end
    usePeak(fmovie(mi)-(j-1))=1;
end
aligned.cueZone_onVoff=usePeak;

figure();
plot(aligned.cueZone,'Color','k');
hold on;
plot(aligned.cue*range(aligned.cueZone)+min(aligned.cueZone),'Color','r');
plot(aligned.cueZone_onVoff*range(aligned.cueZone)+min(aligned.cueZone),'Color','b');
legend({'original from movie','from arduino','fixed'});

end

