function [dprimes,isreaching_out]=plotLearningCurve(exptDir,nameOfCue)

[alltbt,metadata]=combineExptPieces(exptDir,nameOfCue,0.25,1);

[out,alltbt]=getSweepsFromBeh(alltbt);

% Discard any sessions where nth_session is nan
[alltbt,metadata,out]=excludeTrials(alltbt,metadata,out,'nth_session');

isreaching_out=countMouseSuccessDropMiss(alltbt,metadata,out);

% Various methods for calculating dprimes
settings=RTanalysis_settings();
settings.preCueWindow_start=0; % define start of time window from trial onset, in seconds
settings.preCueWindow_end=1.5; % define end of time window from trial onset, in seconds
[dprimes_preCue,hit_rates,FA_rates]=get_dprime_per_session(alltbt,out,metadata,'reachStarts_noPawOnWheel',nameOfCue,settings);

settings=RTanalysis_settings();
settings.preCueWindow_start=3.81; % define start of time window from trial onset, in seconds
settings.preCueWindow_end=5.31; % define end of time window from trial onset, in seconds
[dprimes_postCue,hit_rates,FA_rates]=get_dprime_per_session(alltbt,out,metadata,'reachStarts_noPawOnWheel',nameOfCue,settings);

isreaching_out.dprimes_preCue=dprimes_preCue;
isreaching_out.dprimes_postCue=dprimes_postCue;
dprimes=min([dprimes_preCue; dprimes_postCue],[],1);

end

function out=countMouseSuccessDropMiss(tbt,metadata,ttclass)

% gets the # of times mouse grabbed pellet successfully, dropped pellet or
% missed pellet and as a function of pellet availability

thresh=0.5;
% Get trials where pellet arrived before a cue
%[~,indForPresented]=max(nanmean(tbt.pelletPresented,1)); % when was pellet presented, in indices
[~,indForCue]=max(nanmean(tbt.cueZone_onVoff,1)); % when cue on, indices
pelletPresentOnThisTrial=any(tbt.pelletPresent(:,indForCue-25:indForCue)>thresh,2);

sesstypes=unique(metadata.sessid);
out.hasSuccess=nan(1,length(sesstypes));
out.hasDrop=nan(1,length(sesstypes));
out.hasMiss=nan(1,length(sesstypes));
out.hasReach=nan(1,length(sesstypes));
out.pelletPresent=nan(1,length(sesstypes));
out.touched_pellet=nan(1,length(sesstypes));
for i=1:length(sesstypes)
    currsessid=sesstypes(i);
    out.hasSuccess(i)=sum(any(tbt.reachBatch_success_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
    out.hasDrop(i)=sum(any(tbt.reachBatch_drop_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
    out.hasMiss(i)=sum(any(tbt.reachBatch_miss_reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
    out.hasReach(i)=sum(any(tbt.reachStarts(metadata.sessid==currsessid,indForCue:end)>thresh,2));
    out.pelletPresent(i)=sum(pelletPresentOnThisTrial(metadata.sessid==currsessid));
    out.touched_pellet(i)=sum(ttclass.touched_pellet(metadata.sessid==currsessid));
    out.totalTrials(i)=sum(metadata.sessid==currsessid);
    out.consumedPelletAndNoPawOnWheel(i)=sum(ttclass.consumed_pellet(metadata.sessid==currsessid) & ~ttclass.paw_during_wheel(metadata.sessid==currsessid));
    if isfield(metadata,'optoOnHere')
        out.optoOnHere(i)=nanmean(metadata.optoOnHere(metadata.sessid==currsessid));
    end
    if isfield(metadata,'nth_session')
        out.nth_session(i)=nanmean(metadata.nth_session(metadata.sessid==currsessid));
    end
end

end
