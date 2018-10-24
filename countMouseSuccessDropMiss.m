function out=countMouseSuccessDropMiss(tbt)

% gets the # of times mouse grabbed pellet successfully, dropped pellet or
% missed pellet and as a function of pellet availability

thresh=0.5;
% Get trials where pellet arrived before a cue
%[~,indForPresented]=max(nanmean(tbt.pelletPresented,1)); % when was pellet presented, in indices
[~,indForCue]=max(nanmean(tbt.cueZone_onVoff,1)); % when cue on, indices
pelletPresentOnThisTrial=any(tbt.pelletPresent(:,indForCue-25:indForCue)>thresh,2);

hasSuccess=any(tbt.success_reachStarts(:,indForCue:end)>thresh,2);
hasDrop=any(tbt.drop_reachStarts(:,indForCue:end)>thresh,2);
hasMiss=any(tbt.miss_reachStarts(:,indForCue:end)>thresh,2);
hasReach=any(tbt.reachStarts_pelletPresent(:,indForCue:end)>thresh,2);

out.pelletPresent=pelletPresentOnThisTrial;
out.hasReach=hasReach;
out.hasMiss=hasMiss;
out.hasDrop=hasDrop;
out.hasSuccess=hasSuccess;