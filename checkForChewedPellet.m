function finaldata=checkForChewedPellet(finaldata)

settings=autoReachAnalysisSettings();

minTimePelletChew=settings.chew.minTimeToChewPellet;
withinXSeconds=settings.chew.withinXSeconds;
fps=settings.movie_fps;
% Convert to inds
minIndToPelletChew=floor(minTimePelletChew/(1/fps));
withinXInds=floor(withinXSeconds/(1/fps));

% Check that for each reach classified as a success, there is at least this
% much chewing time after the cue, 

% First for reaches where paw does not start on wheel
finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
[finaldata.success_reachStarts,newDrops]=checkForSufficientChewing(finaldata.success_reachStarts,finaldata.isChewing,minIndToPelletChew,withinXInds);
finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
finaldata.drop_reachStarts(newDrops==1)=1;

% For reaches where paw does start on wheel
finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
[finaldata.success_reachStarts_pawOnWheel,newDrops]=checkForSufficientChewing(finaldata.success_reachStarts_pawOnWheel,finaldata.isChewing,minIndToPelletChew,withinXInds);
finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
finaldata.drop_reachStarts_pawOnWheel(newDrops==1)=1;

end

function [reaches,newDrops]=checkForSufficientChewing(reaches,chewing,minIndToPelletChew,withinXInds)

fi=find(reaches==1);
newDrops=zeros(size(reaches));
for i=1:length(fi)
    currReachInd=fi(i);
    % is there enough chewing within X seconds of this reach
    chewInds=sum(chewing(currReachInd:currReachInd+withinXInds)>0.5);
    if chewInds<minIndToPelletChew % not enough chewing to be consistent with eating pellet
        reaches(currReachInd)=0; % not a successful reach
        newDrops(currReachInd)=1; % actually a drop
    end
end

end
