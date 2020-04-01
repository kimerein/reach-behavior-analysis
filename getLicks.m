function out=getLicks(reachData)

% reachData is intensity in user-defined reach zone
out.rawData=reachData;

% user-set constants
settings=autoReachAnalysisSettings();
maxReachFrequency=settings.reach.maxReachFrequency; % in Hz, the maximum frequency at which mouse can reach
movie_fps=settings.movie_fps; % movie frame rate, in frames per second
reachThresh=settings.reach.reachThresh; % after non-parametric transformation of reachData
holdThreshold=settings.reach.holdThreshold; % in seconds -- if any reach lasts longer than 10 s, this is not a reach -- this is a hold
plotOutput=settings.reach.plotOutput; % if 1, plot output, else do not plot

minReachPeriod=1/maxReachFrequency; % in seconds
movieFramePeriod=1/movie_fps; % in seconds
minIndsBetweenReach=floor(minReachPeriod/movieFramePeriod);

reachData=nonparamZscore(reachData); % non-parametric Z score
isReach=reachData>reachThresh;

% Find peaks
[~,reachPeakLocs,reachWidths]=findpeaks(reachData,'MinPeakDistance',minIndsBetweenReach,'MinPeakHeight',reachThresh,'WidthReference','halfheight');
reachPeaks=zeros(size(reachData));
reachPeaks(reachPeakLocs)=1;
currReachWidth=zeros(size(reachData));
currReachWidth(reachPeakLocs)=reachWidths;

% Find reach stretches lasting longer than holdThreshold
holdIndsThresh=floor(holdThreshold/movieFramePeriod);
isHold=zeros(1,length(reachData));
holdLengths=zeros(1,length(reachData));
currHoldLength=0;
holdStarts=zeros(1,length(reachData));
for i=1:length(isHold)-(holdIndsThresh-1)
    if all(reachData(i:i+(holdIndsThresh-1))>reachThresh)
        isHold(i:i+(holdIndsThresh-1))=1;
        currHoldLength=currHoldLength+1;
    else
        if currHoldLength>0
            holdLengths(i+1-currHoldLength:i+(holdIndsThresh-1))=currHoldLength+holdIndsThresh-1;
            holdStarts(i+1-currHoldLength)=1;
            currHoldLength=0;
        end
    end
end

% If reach is in a hold, delete reach
isReach=single((isReach==1) & (isHold==0));
reachPeaks=single((reachPeaks==1) & (isHold==0));
currReachWidth(isHold==1)=holdLengths(isHold==1);

% Map nans
isReach(isnan(reachData)==1)=nan;
reachPeaks(isnan(reachData)==1)=nan;
currReachWidth(isnan(reachData)==1)=nan;
isHold(isnan(reachData)==1)=nan;
holdStarts(isnan(reachData)==1)=nan;
holdLengths(isnan(reachData)==1)=nan;

reachInds=find(reachPeaks==1);

% Find reach beginnings
reachBegins=nan(1,length(reachInds));
for i=1:length(reachInds)
    % When does this reach begin?
    reachBegins(i)=findReachBeginning(isReach,reachInds(i));
end

% Find reach endings
reachEnds=nan(1,length(reachInds));
for i=1:length(reachInds)
    % When does this reach end?
    reachEnds(i)=findReachEnding(isReach,reachInds(i));
end

out.raw.isReach=isReach;
out.raw.reachPeaks=reachPeaks;
out.raw.currReachWidth=currReachWidth;
out.raw.reachBegins=reachBegins;
out.raw.reachEnds=reachEnds;
out.raw.reachInds=reachInds;

% Consolidate reaches, find first in each stretch
for i=1:length(reachBegins)
    isInStretch=find(reachInds>=reachBegins(i) & reachInds<=reachEnds(i));
    if length(isInStretch)>1
        % Consolidate reach beginnings and endings
        % Set reach ending for this reach to be the ending of the last
        % reach in stretch
        reachEnds(i)=max(reachEnds(isInStretch));
        % More than one reach peak in this stretch -- take first
        % Set later reaches to nan
        reachInds(isInStretch(2:end))=nan;
        reachBegins(isInStretch(2:end))=nan;
        reachEnds(isInStretch(2:end))=nan;
    end
end
reachInds=reachInds(~isnan(reachInds));
reachBegins=reachBegins(~isnan(reachBegins));
reachEnds=reachEnds(~isnan(reachEnds));
reachPeaks(~ismember(1:length(isReach),reachInds) & ~isnan(reachPeaks))=0;
currReachWidth(~ismember(1:length(isReach),reachInds) & ~isnan(currReachWidth))=0;

out.firstReachInds=reachInds;
out.isReach=isReach;
out.reachPeaks=reachPeaks;
out.currReachWidth=currReachWidth;
out.reachBegins=reachBegins;
out.reachEnds=reachEnds;

out.isHold=isHold;
out.holdStarts=holdStarts;
out.holdLengths=holdLengths;

if plotOutput==1
   figure(); 
   plot(reachData);
   hold on;
   plot(out.isReach.*max(reachData),'Color','r');
   plot(out.reachPeaks.*max(reachData),'Color','k');
   plot(out.isHold.*max(reachData),'Color','g');
   leg={'lick zone intensity','is licking','lick peak','is holding'};
   title('Lick Classification');
   legend(leg);
end

end

function out=findReachBeginning(isReach,currReachInd)

if isReach(currReachInd)==0
    % is not reach
    out=nan;
    return
end
for i=currReachInd-1:-1:1
    if isReach(i)==0
        % no longer in this reach
        out=i+1;
        return
    end
end
out=1;

end

function out=findReachEnding(isReach,currReachInd)

if isReach(currReachInd)==0
    % is not reach
    out=nan;
    return
end
for i=currReachInd+1:length(isReach)
    if isReach(i)==0
        % no longer in this reach
        out=i-1;
        return
    end
end
out=length(isReach);

end