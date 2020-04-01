function savehandles=discardLastNFrames(savehandles)

if isfield(savehandles,'discardLastN')
    settings=alignmentSettings(savehandles.discardLastN);
else
    settings=alignmentSettings;
end

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
    elseif ~ischar(savehandles.(f{i}))
        % Event refer to inds
        temp=savehandles.(f{i});
        temp=temp(~(temp>=lastNFrames(1) & temp<=lastNFrames(end)));
        savehandles.(f{i})=temp;
    end
end

savehandles.discardLastN=discardLastN;