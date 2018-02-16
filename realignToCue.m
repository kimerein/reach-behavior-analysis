function realign_tbt=realignToCue(tbt,useAsCue)

% line up cues

% find first cue ind on
cue=tbt.(useAsCue);
[~,cueind]=max(nanmean(cue,1));

f=fieldnames(tbt);
for i=1:length(f)
    realign_tbt.(f{i})=nan(size(tbt.(f{i})));
end

fi=nan(1,size(cue,1));
for i=1:size(cue,1)
    temp=find(cue(i,:)>0.5,1,'first');
    % if cue is missing from this trial, drop this trial
    if isempty(temp)
        fi(i)=nan;
    else
        fi(i)=temp;
    end
end

% realign cue and rest of fields
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
        elseif fi(j)==cueind
            % already aligned
            temp=currfield(j,:);
        elseif fi(j)<cueind
            % shift back in time
            temp=[nan(1,cueind-fi(j)) currfield(j,1:end-(cueind-fi(j)))];
        elseif fi(j)>cueind
            % shift forward in time
            temp=[currfield(j,1+(fi(j)-cueind):end) nan(1,fi(j)-cueind)];
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
            
            
    