function eat=removeLicksFromEat(eat,zoneVals)

eat.isChewing_backup=eat.isChewing;
temp=eat.isChewing;

% When paws or tongue in lick zone, not chewing
% temp=eat.isChewing;
% temp(eat.licks.isReach==1)=0;
% eat.isChewing=temp;

% Only take chewing bouts consistent with chewing A PELLET (i.e., that last
% long enough, longer than settings.chew.minTimeToChewPellet)
settings=autoReachAnalysisSettings();
chewedLongEnough=zeros(size(temp));
binaryVector=temp>0.5; % is chewing
[labeledVector,numRegions]=bwlabel(binaryVector);
% Measure length of each stretch
measurements=regionprops(labeledVector,temp,'Area','PixelIdxList');
for k=1:numRegions
    if measurements(k).Area>=floor(settings.chew.minTimeToChewPellet/(1/settings.movie_fps))
        chewedLongEnough(measurements(k).PixelIdxList)=1;
    end
end
eat.isChewing=chewedLongEnough;
