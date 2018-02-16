function [x,y,x_cue,y_cue]=plotOnlyCuedResponse(tbt,useAsCue,whichReach)

cueshiftback=0.2; % is seconds, time to shift cue onset back, to compensate for resampling aliasing
dofirstreach=0; % if 1, will plot first reach in each trial only (first reach after cue)
excludePawOnWheel=1; % if 1, will exclude reaches when paw begins from wheel

% remove opto trials
if isfield(tbt,'optoOn')
    isOpto=any(tbt.optoOn>0.5,2);
    f=fieldnames(tbt);
    for i=1:length(f)
        temp=tbt.(f{i});
        if size(temp,1)~=length(isOpto)
            % skip
            continue
        end
        tbt.(f{i})=temp(isOpto==0,:);
    end
end
        
% get times
timespertrial=nanmean(tbt.times,1);

% which reaches to plot
temp=tbt.(whichReach);

avcue=nanmean(tbt.(useAsCue),1);
[~,maind]=max(avcue);

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
end

% exclude reaches when paw begins from wheel?
if excludePawOnWheel==1
    temp(tbt.pawOnWheel==1)=0;
end

% plot results
figure();
plot(timespertrial-cueshiftback,nanmean(tbt.(useAsCue),1),'Color','b');
hold on;
plot(timespertrial,smooth(nansum(temp,1)./size(temp,1),7),'Color','k');
x=timespertrial;
y=nansum(temp,1)./size(temp,1);
x_cue=timespertrial-cueshiftback;
y_cue=nanmean(tbt.(useAsCue),1);