function plotCombinedResults(tbt,useAsCue,whichReach)

secondsSolenoidOn=2; % at least this many seconds of solenoid open, per trial, to be considered licking trial
cueshiftback=0.2; % is seconds, time to shift cue onset back, to compensate for resampling aliasing
dofirstreach=1; % if 1, will plot first reach in each trial only (first reach after cue)
excludePawOnWheel=0; % if 1, will exclude reaches when paw begins from wheel
noReachWindow=1; % if mouse does not reach within this time window of cue, consider this a trial where mouse did not reach
discardChewingTrials=0; % if 1, will discard trials when animal was already chewing before cue
discardTooManyReachesTrials=0; % if 1, will discard trials when animal reached more than 

% get times
timespertrial=nanmean(tbt.times,1);
indsSolenoidOn=floor(secondsSolenoidOn/(mode(diff(nanmean(tbt.times,1)))));
indsNoReach=floor(noReachWindow/(mode(diff(nanmean(tbt.times,1)))));

% is licking trial?
isLicking=nansum(tbt.solenoidOn>0.5,2)>indsSolenoidOn;
% plot check
figure();
plot(tbt.solenoidOn(isLicking==1,:)');
figure();
plot(tbt.solenoidOn(isLicking==0,:)');
isLicking=[isLicking(2:end); 0];
disp('this many licking trials');
disp(sum(isLicking==1));
disp('this many reaching trials');
disp(sum(isLicking==0));

% is opto trial?
isOpto=any(tbt.optoOn>0.5,2);
% plot check
figure();
plot(tbt.optoOn(isOpto==1,:)','Color','r');
hold on;
plot(tbt.optoOn(isOpto==0,:)','Color','k');

% animal is grooming during cue?
avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);
% isGrooming=tbt.isGrooming(:,maind)>0.5;
isGrooming=any(tbt.isGrooming(:,maind-50:maind+50)>0.5,2);
disp('this many grooming trials');
disp(sum(isGrooming==1));

% animal reaches too much -- might indicate grooming
if discardTooManyReachesTrials==1
    beginningInds=floor(2/(mode(diff(nanmean(tbt.times,1))))); % first 2 seconds after cue
    temp=tbt.reachStarts;
    tooMuchReaching=zeros(length(isGrooming),1);
    for i=1:size(temp,1)
        fi=find(temp(i,maind:maind+beginningInds)>0.5);
        if length(fi)>4
            tooMuchReaching(i)=1;
        end
    end
    isGrooming=isGrooming | tooMuchReaching;
end

% animal is chewing during cue?
if discardChewingTrials==1
    isChewing=any(tbt.isChewing(:,maind-50:maind-10)>0.5,2);
    disp('this many chewing trials');
    disp(sum(isChewing==1));
    isGrooming=isGrooming | isChewing;
end
    
% which reaches to plot
temp=tbt.(whichReach);

% plot only first reaches?
if dofirstreach==1
    firstreaches=nan(1,size(temp,1));
    % note that everything has already been cue-aligned
    for i=1:size(temp,1)
        fi=find(temp(i,maind+1:end),1,'first')+maind;
        if ~isempty(fi)
            firstreaches(i)=fi;
        end
        temp(i,:)=zeros(size(temp(i,:)));
        temp(i,fi)=1;
    end
    p=ranksum(firstreaches(isLicking==0 & isOpto==0 & isGrooming==0),firstreaches(isLicking==0 & isOpto==1 & isGrooming==0));
    disp('p-value for reaching block');
    disp(p);
    p=ranksum(firstreaches(isLicking==1 & isOpto==0 & isGrooming==0),firstreaches(isLicking==1 & isOpto==1 & isGrooming==0));
    disp('p-value for licking block');
    disp(p);
end

% count misses
fi=nan(size(temp,1));
for i=1:size(temp,1)
    fitemp=find(temp(i,maind+1:end),1,'first')+maind;
    if ~isempty(fitemp)
        fi(i)=fitemp;
    end
end
noreachtrials=fi>maind+indsNoReach;
disp('fraction no reach misses in control reaching');
disp(nansum(noreachtrials(isLicking==0 & isOpto==0 & isGrooming==0))./sum(isLicking==0 & isOpto==0 & isGrooming==0));
disp('fraction no reach misses in opto reaching');
disp(nansum(noreachtrials(isLicking==0 & isOpto==1 & isGrooming==0))./sum(isLicking==0 & isOpto==1 & isGrooming==0));
disp('fraction no reach misses in control licking');
disp(nansum(noreachtrials(isLicking==1 & isOpto==0 & isGrooming==0))./sum(isLicking==1 & isOpto==0 & isGrooming==0));
disp('fraction no reach misses in opto licking');
disp(nansum(noreachtrials(isLicking==1 & isOpto==1 & isGrooming==0))./sum(isLicking==1 & isOpto==1 & isGrooming==0));

% exclude reaches when paw begins from wheel?
if excludePawOnWheel==1
    temp(tbt.pawOnWheel==1)=0;
end

% plot results
doNorm=0;
% timespertrial2=nanmean([timespertrial(1:3:end); timespertrial(2:3:end); timespertrial(3:3:end)],1);
% temp=(temp(:,1:3:end)+temp(:,2:3:end)+temp(:,3:3:end))/3;
% timespertrial2=nanmean([timespertrial(1:2:end); timespertrial(2:2:end)],1);
% temp=(temp(:,1:2:end)+temp(:,2:2:end))/2;
timespertrial2=timespertrial;
% for i=1:size(temp,1)
%     temp(i,:)=smooth(temp(i,:),5);
% end
if doNorm==1
    figure();
    plot(timespertrial-cueshiftback,nanmean(tbt.(useAsCue),1));
    hold on;
    plot(timespertrial2,nansum(temp(isLicking==0 & isOpto==0 & isGrooming==0,:),1)./sum(isLicking==0 & isOpto==0 & isGrooming==0),'Color','c');
    plot(timespertrial2,nansum(temp(isLicking==0 & isOpto==1 & isGrooming==0,:),1)./sum(isLicking==0 & isOpto==1 & isGrooming==0),'Color','m');
    figure(); plot(timespertrial-cueshiftback,nanmean(tbt.cueZone_onVoff,1)); 
    hold on;
    plot(timespertrial2,nansum(temp(isLicking==1 & isOpto==0 & isGrooming==0,:),1)./sum(isLicking==1 & isOpto==0 & isGrooming==0),'Color','k');
    plot(timespertrial2,nansum(temp(isLicking==1 & isOpto==1 & isGrooming==0,:),1)./sum(isLicking==1 & isOpto==1 & isGrooming==0),'Color','r');
    figure(); plot(timespertrial-cueshiftback,nanmean(tbt.cueZone_onVoff,1)); 
    hold on;
    plot(timespertrial2,nansum(temp(isOpto==0 & isGrooming==0,:),1)./sum(isOpto==0 & isGrooming==0),'Color','k');
    plot(timespertrial2,nansum(temp(isOpto==1 & isGrooming==0,:),1)./sum(isOpto==1 & isGrooming==0),'Color','r');
else
    figure();
    plot(timespertrial-cueshiftback,nanmean(tbt.(useAsCue),1));
    hold on;
    plot(timespertrial2,nansum(temp(isLicking==0 & isOpto==0 & isGrooming==0,:),1),'Color','c');
    plot(timespertrial2,nansum(temp(isLicking==0 & isOpto==1 & isGrooming==0,:),1),'Color','m');
    figure(); plot(timespertrial-cueshiftback,nanmean(tbt.cueZone_onVoff,1)); 
    hold on;
    plot(timespertrial2,nansum(temp(isLicking==1 & isOpto==0 & isGrooming==0,:),1),'Color','k');
    plot(timespertrial2,nansum(temp(isLicking==1 & isOpto==1 & isGrooming==0,:),1),'Color','r');
    figure(); plot(timespertrial-cueshiftback,nanmean(tbt.cueZone_onVoff,1)); 
    hold on;
    plot(timespertrial2,nansum(temp(isOpto==0 & isGrooming==0,:),1),'Color','k');
    plot(timespertrial2,nansum(temp(isOpto==1 & isGrooming==0,:),1),'Color','r');
end