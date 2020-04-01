function [tbt,useCues]=discardMisaligned(tbt,useAsCue)

realignToEnd=1; % if 1, will realign to end of cue

cue=tbt.(useAsCue);

avcue=nanmean(cue,1);
[~,maind]=max(avcue);

fi=nan(1,size(cue,1));
for i=1:size(cue,1)
    temp=find(cue(i,maind+1:end)<0.5,1,'first')+maind;
    % if cue is missing from this trial, drop this trial
    if isempty(temp)
        fi(i)=nan;
    else
        fi(i)=temp;
    end
end
mostcommon=mode(fi);

if realignToEnd==1
    useCues=[];
    f=fieldnames(tbt);
    for i=1:length(f)
        currfield=tbt.(f{i});
        newfield=nan(size(currfield));
        if size(currfield,1)~=size(tbt.(useAsCue),1)
            % skip this
            continue
        end
        for j=1:size(cue,1)
            % realign each trial
            if isnan(fi(j))
                % exclude this trial
                % nan out
                temp=nan(size(currfield(j,:)));
            elseif fi(j)==mostcommon
                % already aligned
                temp=currfield(j,:);
            elseif fi(j)<mostcommon
                % shift back in time
                temp=[nan(1,mostcommon-fi(j)) currfield(j,1:end-(mostcommon-fi(j)))];
            elseif fi(j)>mostcommon
                % shift forward in time
                temp=[currfield(j,1+(fi(j)-mostcommon):end) nan(1,fi(j)-mostcommon)];
            end
            newfield(j,:)=temp;
        end
        % save into tbt
        realign_tbt.(f{i})=newfield;
    end
    
    % check alignment
    figure();
    plot(realign_tbt.(useAsCue)');
    title('Re-aligned cues');
    tbt=realign_tbt;
else
    useCues=fi==mostcommon-1;
    
    f=fieldnames(tbt);
    for i=1:length(f)
        if size(tbt.(f{i}),1)~=size(cue,1)
            continue
        end
        temp=tbt.(f{i});
        tbt.(f{i})=temp(useCues,:);
    end
end
