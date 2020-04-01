function finaldata=checkForChewedPellet(finaldata)

settings=autoReachAnalysisSettings();

minTimePelletChew=settings.chew.minTimeToChewPellet;
withinXSeconds=settings.chew.withinXSeconds;
priorSeconds=settings.chew.priorToReach_chewWindow;
dropIfChewingBefore=settings.chew.dropIfChewingBefore;
minTimeMoreStringent=settings.chew.minTimeToChew_afterReach;
fps=settings.movie_fps;
% Convert to inds
minIndToPelletChew=floor(minTimePelletChew/(1/fps));
withinXInds=floor(withinXSeconds/(1/fps));
priorXInds=floor(priorSeconds/(1/fps));
minIndMoreStringent=floor(minTimeMoreStringent/(1/fps));

% Check that for each reach classified as a success, there is at least this
% much chewing time after the cue, 

% First for reaches where paw does not start on wheel
finaldata.success_reachStarts_backup=finaldata.success_reachStarts;
[finaldata.success_reachStarts,newDrops]=checkForSufficientChewing(finaldata.success_reachStarts,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent);
finaldata.drop_reachStarts_backup=finaldata.drop_reachStarts;
finaldata.drop_reachStarts(newDrops==1)=1;

% For reaches where paw does start on wheel
finaldata.success_reachStarts_pawOnWheel_backup=finaldata.success_reachStarts_pawOnWheel;
[finaldata.success_reachStarts_pawOnWheel,newDrops]=checkForSufficientChewing(finaldata.success_reachStarts_pawOnWheel,finaldata.isChewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent);
finaldata.drop_reachStarts_pawOnWheel_backup=finaldata.drop_reachStarts_pawOnWheel;
finaldata.drop_reachStarts_pawOnWheel(newDrops==1)=1;

end

function [reaches,newDrops]=checkForSufficientChewing(reaches,chewing,minIndToPelletChew,withinXInds,dropIfChewingBefore,priorXInds,minIndMoreStringent)

fi=find(reaches==1);
newDrops=zeros(size(reaches));
for i=1:length(fi)
    currReachInd=fi(i);
    % is there enough chewing within X seconds of this reach
    if currReachInd+withinXInds>length(chewing)
        chewInds=sum(chewing(currReachInd:end)>0.5);
    else
        chewInds=sum(chewing(currReachInd:currReachInd+withinXInds)>0.5);
    end
    if chewInds<minIndToPelletChew % not enough chewing to be consistent with eating pellet
        reaches(currReachInd)=0; % not a successful reach
        newDrops(currReachInd)=1; % actually a drop
    end
    if dropIfChewingBefore==1 && newDrops(currReachInd)==0
        % was mouse chewing BEFORE reach?
        if currReachInd-priorXInds<1
            chewInds_before=sum(chewing(1:currReachInd-1)>0.5);
        else
            chewInds_before=sum(chewing(currReachInd-priorXInds:currReachInd-1)>0.5);
        end
        if chewInds_before>floor((minIndToPelletChew/withinXInds)*priorXInds)
            % mouse was chewing before reach
            % did mouse chew long enough after reach, consistent with
            % consumption of full pellet?
            if chewInds<minIndMoreStringent % not enough chewing to be consistent with eating pellet
                reaches(currReachInd)=0; % not a successful reach
                newDrops(currReachInd)=1; % actually a drop
            end
        end
    end
end

end
