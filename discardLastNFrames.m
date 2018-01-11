function savehandles=discardLastNFrames(savehandles)

settings=alignmentSettings;

discardLastN=settings.discardLastN;
if discardLastN==0
    return
end
savehandles.discardLastN=discardLastN;
totalFrames=sum(~isnan(savehandles.LEDvals))+savehandles.discardFirstNFrames;
lastNFrames=totalFrames-discardLastN+1:totalFrames;

f=fieldnames(savehandles);
for i=1:length(f)
    if sum(~isnan(savehandles.(f{i})))==sum(~isnan(savehandles.LEDvals))
        % Each ind is a frame
        lastNotNan=find(~isnan(savehandles.(f{i})),1,'last');
        temp=savehandles.(f{i});
        temp(lastNotNan-discardLastN+1:lastNotNan)=nan;
        savehandles.(f{i})=temp;
    else
        % Event refer to inds
        temp=savehandles.(f{i});
        temp=temp(~(temp>=lastNFrames(1) & temp<=lastNFrames(end)));
        savehandles.(f{i})=temp;
    end
end