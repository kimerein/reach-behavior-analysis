function savehandles=discardLastNFrames(savehandles)

settings=alignmentSettings;

discardLastN=settings.discardLastN;
if discardLastN==0
    return
end
totalFrames=sum(~isnan(savehandles.LEDvals))+savehandles.discardFirstNFrames;
lastNFrames=totalFrames-discardLastN+1:totalFrames;

framesInLED=sum(~isnan(savehandles.LEDvals));

f=fieldnames(savehandles);
for i=1:length(f)
    if sum(~isnan(savehandles.(f{i})))==framesInLED
        % Each ind is a frame
        lastNotNan=find(~isnan(savehandles.(f{i})),1,'last');
        temp=savehandles.(f{i});
        temp(lastNotNan-discardLastN+1:lastNotNan)=nan;
        savehandles.(f{i})=temp;
    elseif ~islogical(savehandles.(f{i})) & ~ischar(savehandles.(f{i}))
        % Event refer to inds
        temp=savehandles.(f{i});
        temp=temp(~(temp>=lastNFrames(1) & temp<=lastNFrames(end)));
        savehandles.(f{i})=temp;
    end
end

% Make other fields match
for i=1:length(f)
    if islogical(savehandles.(f{i}))
        % Refers to reaches
        % Discard values correponding to the discarded reaches
        temp=savehandles.(f{i});
        temp=temp(1:length(savehandles.reachStarts));
        savehandles.(f{i})=temp;
    end
end

savehandles.discardLastN=discardLastN;