function isCorrect=checkBehaviorClassifications(beh_tbt)

% movieframeinds gives index into movie frames
% for a random sampling of different reach type classifications, have user
% check video for correct classification

nTrialsPerType=5; % number of trials of each type to check
% whichTypesToCheck={'reachBatch_success_reachStarts','reachBatch_drop_reachStarts','reachBatch_miss_reachStarts','pelletmissingreach_reachStarts'};
whichTypesToCheck={'reachBatch_drop_reachStarts'};
lowThresh=0.05;

isCorrect=nan(length(whichTypesToCheck),3); % rows are different types of events, first col is number correct, second col is total number checked, third col is number with timing off
for i=1:length(whichTypesToCheck)
    currCheck=whichTypesToCheck{i};
    temp=beh_tbt.(currCheck);
    whichtrials=find(any(temp>lowThresh,2));
    ntotaltrials=length(whichtrials);
    takeN=nTrialsPerType;
    if takeN>ntotaltrials
        takeN=ntotaltrials;
    end
    checkthesetrials=randsample(ntotaltrials,takeN);
    countingcorrect=nan(1,takeN);
    timingoff=nan(1,takeN);
    for j=1:length(checkthesetrials)
        checkingtrial=whichtrials(checkthesetrials(j));
        if ~isfield(beh_tbt,'from_first_video')
           whichvid=0;
        else
            if beh_tbt.from_first_video(checkingtrial,1)>lowThresh
                whichvid=1;
            elseif beh_tbt.from_second_video(checkingtrial,1)>lowThresh
                whichvid=2;
            elseif beh_tbt.from_third_video(checkingtrial,1)>lowThresh
                whichvid=3;
            else
                disp('May be problem -- this trial is not from any video');
                continue
            end
        end
        aretypeframes=find(temp(checkingtrial,:)>lowThresh);
        useit=randsample(length(aretypeframes),1);
        whichframe=beh_tbt.movieframeinds(checkingtrial,aretypeframes(useit));
        opts.Default='No';
        opts.Interpreter='none';
        answer=questdlg(['Movie ' num2str(whichvid) ' frame ' num2str(whichframe) ' is a ' currCheck '?'],...
                        'Checking event classifications',...
                        'Yes','No','Yes but timing off',opts);
        switch answer
            case 'Yes'
                countingcorrect(j)=1;
                timingoff(j)=0;
            case 'Yes but timing off'
                countingcorrect(j)=1;
                timingoff(j)=1;
            case 'No'
                countingcorrect(j)=0;
                timingoff(j)=nan;
        end
        pause;
    end
    isCorrect(i,1)=nansum(countingcorrect);
    isCorrect(i,2)=length(countingcorrect);
    isCorrect(i,3)=nansum(timingoff);
    disp([currCheck ': ' num2str(isCorrect(i,1)) ' of ' num2str(isCorrect(i,2)) ' are correct']);
end

end