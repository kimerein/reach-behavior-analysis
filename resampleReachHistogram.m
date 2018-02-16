function tbt=resampleReachHistogram(tbt,binTimes)

disp('Started resampling for fixed time bins');
% Put tbt into fixed time bins
f=fieldnames(tbt);
for i=1:length(f)
    disp(['Resampling ' f{i}]);
    temp=tbt.(f{i});
    newTimeTbt_temp=nan(size(temp,1),length(binTimes));
    if size(temp,1)~=size(tbt.cue,1)
        continue
    end
    for j=1:size(temp,1)
        temp2=temp(j,:);
        if isnan(temp2(1))
            temp2(~isnan(tbt.times(j,:)))=0;
        end
        % cut off anything that does not match times
        temp2(isnan(tbt.times(j,:)))=nan;
        curr=timeseries(temp2(~isnan(temp2)),tbt.times(j,~isnan(temp2)));
        if all(temp2(~isnan(temp2))==0)
            newTimeTbt_temp(j,:)=zeros(1,length(binTimes));
        else
            res_curr=resample(curr,binTimes);
            newTimeTbt_temp(j,1:length(res_curr.data))=res_curr.data;
        end
    end
    newTimeTbt.(f{i})=newTimeTbt_temp;
end
disp('Finished resampling for fixed time bins');
tbt=newTimeTbt;