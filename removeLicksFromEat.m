function eat=removeLicksFromEat(eat,zoneVals)

% When paws or tongue in lick zone, not chewing
temp=eat.isChewing;
temp(eat.licks.isReach==1)=0;
eat.isChewing=temp;
