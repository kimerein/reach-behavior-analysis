function aligned=getAlignment(out,moviefps,handles,settings)

% Note that microSD (Arduino) output is timed in ms
% Whereas video is timed in frames per sec

% Get settings for alignment
if isempty(settings)
    settings=alignmentSettings();
end

% Remove incomplete reach detections
if length(handles.eatTime)<length(handles.reachStarts)
    handles.eatTime=[handles.eatTime ones(1,length(handles.reachStarts)-length(handles.eatTime))*handles.eatTime(end)];
end
isnotnaninds=~isnan(mean([handles.reachStarts' handles.pelletTime' handles.eatTime' handles.pelletMissing'],2));
f=fieldnames(handles);
l=length(handles.reachStarts);
for i=1:length(f)
    if length(handles.(f{i}))==l
        temp=handles.(f{i});
        handles.(f{i})=temp(isnotnaninds);
    end
end

% Try to align based on distractor LED from movie and Arduino output
temp_LED=handles.LEDvals;
temp_LED=temp_LED(~isnan(temp_LED));
% temp_LED=temp_LED-smooth(temp_LED,51)';
% temp_LED=temp_LED-min(temp_LED);
% [~,whtorm]=rmoutliers(temp_LED,2,"ThresholdFactor",20);
whtorm=temp_LED>20*mad(temp_LED,[],2)+median(temp_LED,2) | temp_LED<median(temp_LED,2)-20*mad(temp_LED,[],2);
temp_LED(whtorm)=prctile(temp_LED,10);
threshForOnVsOff=min(temp_LED)+settings.fractionRange*range(temp_LED);
% threshForOnVsOff=nanmean([max(temp_LED) min(temp_LED)])-0.4*(max(temp_LED)-min(temp_LED));

% Check for single points above thresh -- LED duration is more than 1 frame
isGreater=temp_LED>threshForOnVsOff;
justOnePoint=[diff(isGreater(1:2))==-1 diff(isGreater(1:end-1))==1 & diff(isGreater(2:end))==-1];
temp_LED(justOnePoint)=threshForOnVsOff-1;

figz(1)=figure();
movie_times=0:(1/moviefps)*1000:(length(temp_LED)-1)*((1/moviefps)*1000);
plot(movie_times,temp_LED,'Color','b');
hold on;
line([0 (length(temp_LED)-1)*((1/moviefps)*1000)],[threshForOnVsOff threshForOnVsOff],'Color','r');
title('Threshold for distinguishing LED on vs off');

% Get when LED was on in movie vs off
movie_LED=single(temp_LED>threshForOnVsOff);
formovieframeinds_align=movie_LED;
movie_LED_for_finalalignment=movie_LED;

% Find best alignment of distractor LED in movie and Arduino output -- note
% different sampling rates
temp=out.distractorOn';
testRunLED=out.distractorOn;

temptimes=[];
for i=1:size(out.allTrialTimes,1)
    temptimes=[temptimes out.allTrialTimes(i,:)];
end
temp=temp(1:end);
arduino_LED=temp(~isnan(temptimes));
arduino_times=temptimes(~isnan(temptimes));
backup_arduino_times=arduino_times;

% Find alignment
% First down-sample arduino LED
arduino_dec=settings.arduino_dec;
movie_dec=settings.movie_dec;

arduino_LED=decimate(arduino_LED,arduino_dec);
arduino_times=decimate(arduino_times,arduino_dec);

testRun_movieLED=double(movie_LED);

movie_LED=decimate(double(movie_LED),movie_dec);
movie_times=decimate(movie_times,movie_dec);

% Do an initial alignment
temp=arduino_LED;
arduino_LED(temp>=0.5)=1;
arduino_LED(temp<0.5)=0;
temp=movie_LED;
movie_LED(temp>=0.5)=1;
movie_LED(temp<0.5)=0;

% Discard beginning of Arduino LED
arduino_LED(arduino_times<settings.discardTimeArduinoLED*1000)=nan; % convert discardTimeArduinoLED to ms

% Throw out LED distractor on intervals less than settings.useDistractorThresh
% This deals with skipping of low frame rate DVR
allEvents_arduino_LED=arduino_LED;
allEvents_movie_LED=movie_LED;
arduino_LED=throwOutOnStretches(arduino_LED,arduino_times);
[movie_LED,throwOutMovie]=throwOutOnStretches(movie_LED,movie_times);

[pks_arduino,locs_arduino]=findpeaks(arduino_LED);
arduino_LED_ITIs=diff(arduino_times(locs_arduino));
[pks,locs]=findpeaks(movie_LED);
movie_LED_ITIs=diff(movie_times(locs)); 

if settings.isInSecondHalf==true
    backup_arduino_LED_ITIs=arduino_LED_ITIs;
    midLength=floor(settings.fractionThroughArduino*length(arduino_LED_ITIs))+20;
    arduino_LED_ITIs=arduino_LED_ITIs(midLength+1:end);
end
temp1=arduino_LED_ITIs./max(arduino_LED_ITIs);
temp2=movie_LED_ITIs./max(movie_LED_ITIs);
% temp1=temp1-nanmean(temp1);
% temp2=temp2-nanmean(temp2);
if isempty(settings.maxlagForInitialAlign) 
    [X,Y,D]=alignsignals(temp1,temp2); 
else
    [X,Y,D]=alignsignals(temp1,temp2,settings.maxlagForInitialAlign); 
end
if settings.isInSecondHalf==true
    D=D-midLength;
    arduino_LED_ITIs=backup_arduino_LED_ITIs;
    if D>0
        X=[zeros(1,D) arduino_LED_ITIs./max(arduino_LED_ITIs)];
        Y=movie_LED_ITIs./max(movie_LED_ITIs);
    elseif D<0
        Y=[zeros(1,-D) movie_LED_ITIs./max(movie_LED_ITIs)];
        X=arduino_LED_ITIs./max(arduino_LED_ITIs);
    end
end

tryinc=settings.tryinc; % this is the increment for trying different scalings of movie onto arduino data
if D>0
    error('Why does movie start before Arduino?');
else
    movie_peakloc=1;
    arduino_peakloc=abs(D)+1;
    movie_peak_indexIntoMovie=locs(movie_peakloc);
    arduino_peak_indexIntoArduino=locs_arduino(arduino_peakloc);
    if (-D)+1+length(movie_LED_ITIs)>length(locs_arduino)
        size_of_arduino=length(arduino_LED(locs_arduino((-D)+1):end));
    else
        size_of_arduino=length(arduino_LED(locs_arduino((-D)+1):locs_arduino((-D)+1+length(movie_LED_ITIs))));
    end
    size_of_movie=length(movie_LED(locs(1):locs(end)));
%     size_of_movie=length(movie_LED(locs(383-312):locs(403-312)));
%     size_of_arduino=length(arduino_LED(locs_arduino(383):locs_arduino(403)));
%     movie_peak_indexIntoMovie=locs(383-312);
%     arduino_peak_indexIntoArduino=locs_arduino(383);
    guess_best_scale=size_of_arduino/size_of_movie;
    guess_best_scale1=guess_best_scale;
    % Adjust according to guess_best_scale
    movledinds=1:length(movie_LED);
    movie_LED=resample(movie_LED,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
    movledinds=resample(movledinds,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
    [~,mi]=min(abs(movledinds-movie_peak_indexIntoMovie));
    movie_peak_indexIntoMovie=mi;
    guess_best_delay=arduino_peak_indexIntoArduino-movie_peak_indexIntoMovie;
    trydelays=guess_best_delay+settings.try_delay1:guess_best_delay+settings.try_delay2;
    % Note that fixed, so now best scale is 1
    guess_best_scale=1;
    tryscales=guess_best_scale+settings.try_scale1:tryinc:guess_best_scale+settings.try_scale2;
end

figz(2)=figure();
plot(X,'Color','b');
hold on;
plot(Y,'Color','r');
title('Preliminary alignment of movie distractor intervals onto arduino distractor intervals');
legend({'Arduino distractor intervals','Movie distractor intervals'});

figz(3)=figure();
plot(arduino_LED,'Color','b');
hold on;
plot([nan(1,guess_best_delay) movie_LED],'Color','r');
title('Preliminary alignment of movie distractor onto arduino distractor');
legend({'Arduino distractor','Movie distractor'});

% Wait for user to confirm preliminary alignment
if settings.isOrchestra~=1
    pause;
end
 
% Test signal alignment and scaling
disp('Now refining alignment ...');
sumdiffs=nan(length(tryscales),length(trydelays));
if settings.alignWithAllEvents==1
    backup_movie_LED=allEvents_movie_LED;
    forSecondaryAlignment=backup_movie_LED;
    backup_arduino_LED=allEvents_arduino_LED;
    % Remove LED distractor on intervals that are too short (i.e., may have
    % been missed in movie, shorter than 3 movie frames)
    [backup_movie_LED,throwOutMovie]=throwOutOnStretches(backup_movie_LED,1:length(backup_movie_LED),3);
    backup_movie_LED=resample(backup_movie_LED,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
    throwOutMovie=interp(throwOutMovie,movie_dec);
    movie_LED_for_finalalignment(throwOutMovie>0.5)=0;
    figz(4)=figure(); 
    plot(allEvents_movie_LED,'Color','b'); 
    hold on; 
    plot(movie_LED_for_finalalignment,'Color','r');
    legend({'before removal','after removal'});
    title('Remove short, skipped frame LED distractors before alignment');
else
    figz(4)=figure;
    backup_movie_LED=movie_LED;
    forSecondaryAlignment=backup_movie_LED;
    backup_arduino_LED=arduino_LED;
end
for j=1:length(tryscales)
    if mod(j,10)==0
        disp('Processing ...');
%         disp(j);
    end
    currscale=tryscales(j);
    movie_LED=resample(backup_movie_LED,floor(currscale*(1/tryinc)),floor((1/tryinc)));
    for i=1:length(trydelays)
        currdelay=trydelays(i);
        if mod(i,500)==0
%             disp(i);
        end 
        temp_movie=[nan(1,currdelay) movie_LED];
        temp_arduino=[arduino_LED nan(1,length(temp_movie)-length(arduino_LED))];
        if length(temp_arduino)>length(temp_movie)
            temp_movie=[temp_movie nan(1,length(temp_arduino)-length(temp_movie))];
            sumdiffs(j,i)=nansum(abs(temp_movie-temp_arduino));
        elseif length(temp_movie)>length(temp_arduino)
            % This should not happen
            % Arduino should include movie
            sumdiffs(j,i)=inf;
            error('Why is movie longer than Arduino?');
        else
            sumdiffs(j,i)=nansum(abs(temp_movie-temp_arduino)); 
        end
    end
end
sumdiffs(isinf(sumdiffs))=nan;
ma=max(sumdiffs(:));
sumdiffs(isnan(sumdiffs))=3*ma;
[minval,mi]=min(sumdiffs(:));
[mi_row,mi_col]=ind2sub(size(sumdiffs),mi);

figz(5)=figure(); 
imagesc(sumdiffs);
title('Finding best alignment');
xlabel('Trying different delays');
ylabel('Trying different scales');
disp('Best delay column');
disp(mi_col);
disp('Best scale row');
disp(mi_row);

temp=resample(forSecondaryAlignment,floor(mod(size_of_arduino/size_of_movie,1)*100)+floor((guess_best_scale*100)/100)*100,100);
frontShift=trydelays(mi_col);
scaleBy=tryscales(mi_row);
resampFac=1/tryinc;
best_movie=[nan(1,trydelays(mi_col)) resample(temp,floor(tryscales(mi_row)*(1/tryinc)),floor(1/tryinc))];
% best_movie=[nan(1,trydelays(mi_col)) resample(backup_movie_LED,floor(tryscales(mi_row)*(1/tryinc)),floor(1/tryinc))];
shouldBeLength=length(best_movie);
best_arduino=[backup_arduino_LED nan(1,length(best_movie)-length(backup_arduino_LED))];
movieToLength=length(best_arduino);
if length(best_arduino)>length(best_movie)
    best_movie=[best_movie nan(1,length(best_arduino)-length(best_movie))];
end
figz(6)=figure();
plot(best_movie,'Color','r');
hold on;
plot(best_arduino,'Color','b');
title('Alignment of movie distractor onto arduino distractor before small segment adjustments');
legend({'Movie','Arduino'});

% Allow user to zero out any spurious LED distractors before further
% alignment
if settings.doManualZeroOut==1
    [best_movie,best_arduino]=zeroOutDialog(best_movie,best_arduino);
end

% Then re-align sub-sections of movie to arduino code
alignSegments=settings.alignSegments; % in number of indices
mov_distractor=[];
arduino_distractor=[];
firstInd=find(~isnan(best_movie) & ~isnan(best_arduino),1,'first');
lastBoth=min([find(~isnan(best_movie),1,'last') find(~isnan(best_arduino),1,'last')]);
if lastBoth-firstInd<alignSegments
    alignSegments=ceil((lastBoth-firstInd)/2);
end
segmentInds=firstInd:floor(alignSegments/2):lastBoth;
allMatchedInds=firstInd:lastBoth;
mov_distractor=[mov_distractor nan(1,firstInd-1)];
arduino_distractor=[arduino_distractor nan(1,firstInd-1)];
segmentDelays=nan(1,length(segmentInds));
addZeros_movie=nan(1,length(segmentInds));
addZeros_arduino=nan(1,length(segmentInds));
moveChunks=nan(length(segmentInds),2);
tookTheseIndsOfTemp1=floor(0.25*alignSegments):ceil(0.75*alignSegments);
haveDoneTheseInds=zeros(size(firstInd:lastBoth));
backup_tookTheseIndsOfTemp1=tookTheseIndsOfTemp1;
for i=1:length(segmentInds)-1
    currInd=segmentInds(i);
    if currInd+alignSegments-1>length(best_movie) || currInd+alignSegments-1>length(best_arduino)
        [temp1,temp2,D]=alignsignals(best_movie(currInd:end),best_arduino(currInd:end));
    else
        [temp1,temp2,D]=alignsignals(best_movie(currInd:currInd+alignSegments-1),best_arduino(currInd:currInd+alignSegments-1));
    end
    segmentDelays(i)=D;
    if length(temp1)>length(temp2)
        addZeros_arduino(i)=length(temp1)-length(temp2);
        temp2=[temp2 temp2(end)*ones(1,length(temp1)-length(temp2))];
        addZeros_movie(i)=0;
    elseif length(temp2)>length(temp1)
        addZeros_movie(i)=length(temp2)-length(temp1);
        temp1=[temp1 temp1(end)*ones(1,length(temp2)-length(temp1))];
        addZeros_arduino(i)=0;
    else
        addZeros_movie(i)=0;
        addZeros_arduino(i)=0;
    end
    % Take middle of aligned segments, because alignment at middle tends to
    % be better than alignment at edges
    temp_startInd=find(haveDoneTheseInds==0,1,'first');
    currIndices=currInd:(currInd+alignSegments-1);
    tookTheseIndsOfTemp1(1)=allMatchedInds(temp_startInd)-currIndices(1)+1;
    tookTheseIndsOfTemp1(end)=backup_tookTheseIndsOfTemp1(end);
    temp_endInd=temp_startInd+(tookTheseIndsOfTemp1(end)-tookTheseIndsOfTemp1(1));
    if i==1
        haveDoneTheseInds(1:temp_endInd)=1;
    else
        haveDoneTheseInds(temp_startInd:temp_endInd)=1;
    end
    if i==1
        startAt=1;
        if D>0 % temp1 has been delayed by D samples
            endAt=tookTheseIndsOfTemp1(end)+D;
        else % temp2 has been delayed by D samples
            endAt=tookTheseIndsOfTemp1(end);
        end
    elseif i==length(segmentInds)-1
        % alignment easily messed up at end -- just use delay from previous
        % segment
        if segmentDelays(i-1)>0
            temp1=[ones(1,segmentDelays(i-1))*temp1(1) temp1];
        else
            temp2=[ones(1,-segmentDelays(i-1))*temp2(1) temp2];
        end
        segmentDelays(i)=segmentDelays(i-1);
        D=segmentDelays(i);
        if D>0 % temp1 has been delayed by D samples
            startAt=tookTheseIndsOfTemp1(1)+D;
        else % temp2 has been delayed by D samples
            startAt=tookTheseIndsOfTemp1(1);
        end
        endAt=min([length(temp1) length(temp2)]);
    else
        if D>0 % temp1 has been delayed by D samples
            startAt=tookTheseIndsOfTemp1(1)+D;
            endAt=tookTheseIndsOfTemp1(end)+D;
        else % temp2 has been delayed by D samples
            startAt=tookTheseIndsOfTemp1(1);
            endAt=tookTheseIndsOfTemp1(end);
        end
    end
    % Check to be sure that edges are LED off
    tryj=1;
    while temp1(startAt)>0.5 || temp2(startAt)>0.5
        if temp1(startAt+tryj)<0.5 && temp2(startAt+tryj)<0.5
            startAt=startAt+tryj;
            break
        end
        tryj=tryj+1;
    end
    tryj=1;
    if endAt>length(temp1) || endAt>length(temp2)
        endAt=min([length(temp1) length(temp2)]);
    end
    while temp1(endAt)>0.5 || temp2(endAt)>0.5
        if temp1(endAt-tryj)<0.5 && temp2(endAt-tryj)<0.5
            endAt=endAt-tryj;
            break
        end
        tryj=tryj+1;
    end
    temp1=temp1(startAt:endAt);
    temp2=temp2(startAt:endAt);
    moveChunks(i,1)=startAt;
    moveChunks(i,2)=endAt;
    mov_distractor=[mov_distractor temp1];
    arduino_distractor=[arduino_distractor temp2];
end
figz(7)=figure();
plot(mov_distractor,'Color','r');
hold on;
plot(arduino_distractor,'Color','b');
title('Alignment of movie distractor onto arduino distractor after small segment adjustments');
legend({'Movie','Arduino'});
aligned.movie_distractor=mov_distractor;
aligned.arduino_distractor=arduino_distractor;

% Times from arduino
isTimes=true;
timesfromarduino=alignLikeDistractor(double(backup_arduino_times),0,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks,isTimes); 
aligned.timesfromarduino=timesfromarduino; % note that this is in ms

% Get rid of extraneous arduino LED distractor on intervals, in case of
% skipping these in movie due to low DVR frame rate with respect to LED on
% duration
fi=find(~isnan(mov_distractor),1,'first');
minLEDinterval=settings.minLEDinterval*1000; % convert from seconds to ms
indInterval=nanmean(diff(aligned.timesfromarduino));
indThresh=floor(minLEDinterval/indInterval);
for i=fi:length(arduino_distractor)
    if arduino_distractor(i)>0.5
        % arduino distractor on
        % check for movie distractor on within indThresh of fi
        if i-indThresh<1
            tryInds=1:i+indThresh;
        elseif i+indThresh>length(mov_distractor) || i+indThresh>length(arduino_distractor)
            tryInds=i-indThresh:min([length(mov_distractor) length(arduino_distractor)]);
        else
            tryInds=i-indThresh:i+indThresh;
        end
        if any(mov_distractor(tryInds)>0.5)
            % movie picked up this distractor
        else
            % movie failed to pick up this distractor
            % zero out this on interval in arduino distractor
            arduino_distractor(tryInds)=0;
        end
    end
end
aligned.arduino_distractor=arduino_distractor;
figz(8)=figure();
plot(mov_distractor,'Color','r');
hold on;
plot(arduino_distractor,'Color','b');
title('Alignment of movie distractor onto arduino distractor after removing skipped arduino on intervals');
legend({'Movie','Arduino'});

% Align other signals in same fashion as LED distractor
% From Arduino

% Make cue ONLY the *start* of cue
for i=1:size(out.cueOn,1)
    temp=zeros(size(out.cueOn(i,:)));
    temp(isnan(out.cueOn(i,:)))=nan;
    temp(find(out.cueOn(i,:)>0.5,1,'first'))=1;
    out.cueOn(i,:)=temp;
end

if isfield(out,'falseCueOn')
    for i=1:size(out.falseCueOn,1)
        temp=zeros(size(out.falseCueOn(i,:)));
        temp(isnan(out.falseCueOn(i,:)))=nan;
        temp(find(out.falseCueOn(i,:)>0.5,1,'first'))=1;
        out.falseCueOn(i,:)=temp;
    end
end

% Align cue and/or false cue
aligned.cue=alignDataLikeDistractor(out.cueOn,temptimes,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks);
if isfield(out,'falseCueOn')
    aligned.falseCueOn=alignDataLikeDistractor(out.falseCueOn,temptimes,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks);
end
% Align distractor from arduino
% This is just a test for alignDataLikeDistractor method
aligned.testRunDistractor=alignDataLikeDistractor(testRunLED,temptimes,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks);

aSet=arduinoSettings();
if aSet.noInterlock==1
      out.interlockSamples=zeros(size(out.interlockSamples));
end

% Align other fields from arduino or movie
for i=1:length(settings.alignField)
    n=settings.alignField(i).name;
    if settings.alignField(i).fromarduino==1
        % align like arduino
        aligned.(n)=alignDataLikeDistractor(out.(n),temptimes,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks);
    else
        % align like movie
        aligned.(n)=alignDataLikeMovie(handles.(n),movie_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_movie,scaleBy,resampFac,moveChunks,guess_best_scale1,false);
    end
end

% From movie

% Movie frame inds
temp=1:length(handles.LEDvals);
movieframeinds_raw=double(handles.discardFirstNFrames+temp);
isTimes=true;
movieframeinds=alignDataLikeMovie(movieframeinds_raw,movie_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_movie,scaleBy,resampFac,moveChunks,guess_best_scale1,isTimes);
% movieframeinds=alignLikeDistractor(movieframeinds_raw,1,movie_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_movie,scaleBy,resampFac,moveChunks,guess_best_scale1);
movieframeinds_backup=movieframeinds;

% Align distractor from movie
% This is just a test for alignDataLikeMovie method
aligned.testRunDistractor_movie=alignDataLikeMovie(testRun_movieLED,movie_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_movie,scaleBy,resampFac,moveChunks,guess_best_scale1,false);

% In case looking at second half of savehandles
f=fieldnames(aligned);
temp=aligned.movie_distractor;
for i=1:length(f)
    curr=aligned.(f{i});
    if abs(length(curr)-length(temp))<50
        % fix completely a hack, not good KR
        if length(curr)<length(temp)
            curr=[curr nan(1,abs(length(curr)-length(temp)))];
        end
    end
    aligned.(f{i})=curr(~isnan(temp));
end
if abs(length(movieframeinds)-length(temp))<50
    % fix completely a hack, not good KR
    if length(movieframeinds)<length(temp)
        movieframeinds=[movieframeinds nan(1,abs(length(movieframeinds)-length(temp)))];
    end
end
movieframeinds=movieframeinds(~isnan(temp));
  
% Re-align movie frame inds based on alignment to non-interped LED vals

% Throw out LED distractor on intervals less than settings.useDistractorThresh
% This is one way to deal with skipping of low frame rate DVR
% handles.LEDvals=movie_LED_for_finalalignment;

maxNFramesForLEDtoChange=settings.maxNFramesForLEDtoChange;
deriv_LEDvals=[diff(double(formovieframeinds_align)) 0];
% deriv_LEDvals=[diff(handles.LEDvals) 0];
deriv_thresh=0.5*range(deriv_LEDvals(deriv_LEDvals>=0));
% deriv_thresh=(max(handles.LEDvals)/2)/(maxNFramesForLEDtoChange+1);
[~,peakLocs]=findpeaks(deriv_LEDvals,'MinPeakHeight',deriv_thresh,'MinPeakDistance',3*maxNFramesForLEDtoChange);
% [pks,locs]=findpeaks(deriv_LEDvals);
% peakLocs=locs(pks>deriv_thresh);
% tooclose=diff(peakLocs);
% peakLocs=peakLocs(tooclose>=3*maxNFramesForLEDtoChange);
[~,troughLocs]=findpeaks(-deriv_LEDvals,'MinPeakHeight',deriv_thresh,'MinPeakDistance',3*maxNFramesForLEDtoChange);
% [pks,locs]=findpeaks(-deriv_LEDvals);
% troughLocs=locs(pks>deriv_thresh);
% tooclose=diff(troughLocs);
% troughLocs=troughLocs(tooclose>=3*maxNFramesForLEDtoChange);
rawmovieinds_onto_rescaled=nan(size(movieframeinds));
if movieframeinds_raw(peakLocs(1))<movieframeinds_raw(troughLocs(1))
    % LED first increases
else
    % LED first decreases
    % Throw out first decrease
    troughLocs=troughLocs(2:end);
end
if length(peakLocs)>length(troughLocs)
    peakLocs=peakLocs(1:length(troughLocs));
end
if movieframeinds_raw(peakLocs(1))>=movieframeinds_raw(troughLocs(1))
    error('Problem in raw LED values -- should always turn on, then off');
end
rescaled_thresh=0.5;
temp=aligned.movie_distractor;
[~,f]=findpeaks(temp,'MinPeakHeight',0.08,'MinPeakDistance',indThresh);
temp(f)=1;
deriv_LEDvals_rescaled=[diff(double(temp>0.5)) 0];
% rescaled_thresh=0.25*range(deriv_LEDvals_rescaled(deriv_LEDvals_rescaled>=0));
[~,peakLocs_rescaled]=findpeaks(deriv_LEDvals_rescaled,'MinPeakHeight',rescaled_thresh,'MinPeakDistance',3*maxNFramesForLEDtoChange);
% [pks,locs]=findpeaks(deriv_LEDvals_rescaled);
% peakLocs_rescaled=locs(pks>rescaled_thresh);
% tooclose=diff(peakLocs_rescaled);
% peakLocs_rescaled=peakLocs_rescaled(tooclose>=3*maxNFramesForLEDtoChange);
% peakTimes_rescaled=movieframeinds(peakLocs_rescaled);
[~,troughLocs_rescaled]=findpeaks(-deriv_LEDvals_rescaled,'MinPeakHeight',rescaled_thresh,'MinPeakDistance',3*maxNFramesForLEDtoChange);
% [pks,locs]=findpeaks(-deriv_LEDvals_rescaled);
% troughLocs_rescaled=locs(pks>rescaled_thresh);
% tooclose=diff(troughLocs_rescaled);
% troughLocs_rescaled=troughLocs_rescaled(tooclose>=3*maxNFramesForLEDtoChange);
% troughTimes_rescaled=movieframeinds(troughLocs_rescaled);
% if peakTimes_rescaled(1)<troughTimes_rescaled(1)
if peakLocs_rescaled(1)<troughLocs_rescaled(1)
    % LED first increases
else
    % LED first decreases
    % Throw out first decrease
    troughLocs_rescaled=troughLocs_rescaled(2:end);
%     troughTimes_rescaled=troughTimes_rescaled(2:end);
end
if length(peakLocs_rescaled)>length(troughLocs_rescaled)
    peakLocs_rescaled=peakLocs_rescaled(1:length(troughLocs_rescaled));
%     peakTimes_rescaled=peakTimes_rescaled(1:length(troughLocs_rescaled));
end
k=1;
donotdoalign=settings.donotdoalign;
if length(peakLocs_rescaled)~=length(peakLocs)
    fortitle='May be a problem: different numbers of LED distractor flashes during alignment of movieframeinds';
    disp(fortitle);
else
    fortitle='Good: matching LED distractor flashes during alignment of movieframeinds';
    disp(fortitle);
end
for i=1:length(peakLocs)
    up=movieframeinds_raw(peakLocs(i));
    down=movieframeinds_raw(troughLocs(i));
%     [~,mi_up]=min(abs(peakTimes_rescaled-up));
%     [~,mi_down]=min(abs(troughTimes_rescaled-down));
%     if (k>length(troughLocs_rescaled)) || (k>length(peakLocs_rescaled))
    if (k>length(troughLocs_rescaled)) 
        break
%         donotdoalign=1;
%         disp('donotdoalign');
%         break
    end
    rawmovieinds_onto_rescaled(peakLocs_rescaled(k):troughLocs_rescaled(k))=linspace(up,down,troughLocs_rescaled(k)-peakLocs_rescaled(k)+1);
    k=k+1;
end
rawmovieinds_onto_rescaled=rawmovieinds_onto_rescaled+0.5; % accounts for shift from [diff(handles.LEDvals) 0];
if donotdoalign==0
    % Fill in nans accordingly
    firstnan=find(isnan(rawmovieinds_onto_rescaled),1,'first');
    temp=rawmovieinds_onto_rescaled;
    temp(1:firstnan)=nan;
    firstnotnan=find(~isnan(temp),1,'first');
    safety_counter=1;
    rawframestart=movieframeinds_raw(1);
    while ~isempty(firstnan)
        if safety_counter>20*10^4
            break
        end
        if isempty(rawmovieinds_onto_rescaled(firstnotnan))
            temp=linspace(rawframestart,movieframeinds_raw(end),length(rawmovieinds_onto_rescaled)-firstnan+2);
            rawmovieinds_onto_rescaled(firstnan:length(rawmovieinds_onto_rescaled))=temp(2:end);
            break
        end
        temp=linspace(rawframestart,rawmovieinds_onto_rescaled(firstnotnan),firstnotnan-firstnan+2);
        rawmovieinds_onto_rescaled(firstnan:firstnotnan-1)=temp(2:end-1);
        % Increment
        firstnan=find(isnan(rawmovieinds_onto_rescaled),1,'first');
        rawframestart=rawmovieinds_onto_rescaled(firstnan-1);
        temp=rawmovieinds_onto_rescaled;
        temp(1:firstnan)=nan;
        firstnotnan=find(~isnan(temp),1,'first');
        safety_counter=safety_counter+1;
    end
    aligned.movieframeinds=rawmovieinds_onto_rescaled;
else
    aligned.movieframeinds=movieframeinds;
%     aligned.movieframeinds=smooth(movieframeinds,100);
end

% Plot results
figz(9)=figure();
ha=tight_subplot(7,1,[0.06 0.03],[0.08 0.1],[0.1 0.01]);
currha=ha(1);
axes(currha);
plot(aligned.cue,'Color','r');
xlabel('Cue');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');
title('Results of alignment');

currha=ha(2);
axes(currha);
plot(aligned.movie_distractor,'Color','b');
hold on;
plot(aligned.arduino_distractor,'Color','r');
xlabel('Distractor');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');

currha=ha(3);
axes(currha);
plot(aligned.movieframeinds,'Color','b');
xlabel('Movie frames');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');

currha=ha(4);
axes(currha);
plot(aligned.pelletLoaded,'Color','r');
xlabel('Pellet loaded');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');

currha=ha(5);
axes(currha);
plot(aligned.pelletPresented,'Color','r');
xlabel('Pellet presented');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');

currha=ha(6);
axes(currha);
plot(aligned.timesfromarduino./1000,'Color','r');
xlabel('Times from arduino');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');

currha=ha(7);
axes(currha);
plot(aligned.movieframeinds.*(1/moviefps),'Color','b');
xlabel('Times from movie');
set(currha,'XTickLabel','');
% set(currha,'YTickLabel','');

% Plot realigned movieframeinds
figz(10)=figure();
plot(aligned.movieframeinds);
title(fortitle);

if settings.isOrchestra==1
    endofVfname=regexp(handles.filename,'\.');
    savefig(figz,[handles.filename(1:endofVfname(end)-1) '_alignmentFigs.fig'],'compact');
    close all
end

end

function [movie,arduino]=zeroOutDialog(movie,arduino)

disp('Enter any ranges to zero out in arduino. Enter as vectors in a cell array (e.g., {1:2, 3:4}).');
s=input('Arduino zero out. Enter here: ');
for i=1:length(s)
   arduino(s{i})=0;
end

disp('Enter any ranges to zero out in movie. Enter as vectors in a cell array (e.g., {1:2, 3:4}).');
s=input('Movie zero out. Enter here: ');
for i=1:length(s)
   movie(s{i})=0;
end

end

function [out,throwOut]=throwOutOnStretches(varargin)

data=varargin{1};
times=varargin{2};
if length(varargin)>2
    indThresh=varargin{3};
end

% Get thresh for on vs off
threshForData=nanmin(data)+0.25*range(data);

timeInt=mode(diff(times));
settings=alignmentSettings();
if length(varargin)>2
    threshInds=indThresh;
else
    thresh=settings.useDistractorThresh; % in ms
    threshInds=floor(thresh/timeInt);
end
temp=data;
indon=temp>threshForData;
throwOut=zeros(size(data));
safeguard=1;
while any(indon) && safeguard<length(data)
    fi=find(temp>threshForData,1,'first'); % find first on
    fi2=find(temp(fi:end)<threshForData,1,'first'); % find next off
    if fi2<threshInds
        throwOut(fi:fi+fi2-1)=1;
    end
    temp(1:fi+fi2-1)=0;
    indon=temp>threshForData;
    safeguard=safeguard+1;
end
out=data;
out(throwOut>0.5)=nanmin(data);

end

function X = naninterp(X) 

% Interpolate over NaNs 
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)), 'cubic'); 
return 

end

function out=alignDataLikeMovie(data,movie_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_movie,scaleBy,resampFac,moveChunks,guess_best_scale1,isTimes)

out=alignLikeDistractor(data,1,movie_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_movie,scaleBy,resampFac,moveChunks,guess_best_scale1,isTimes);

end

function out=alignDataLikeDistractor(data,times,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks)

temp=data';
temp=temp(1:end);
testRunDistractor=temp(~isnan(times));
testRunDistractor=alignLikeDistractor(testRunDistractor,0,arduino_dec,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros_arduino,scaleBy,resampFac,moveChunks,[],false); 
out=testRunDistractor./nanmax(testRunDistractor);

end

function outsignal=alignLikeDistractor(signal,scaleThisSignal,decind,frontShift,shouldBeLength,movieToLength,alignSegments,segmentInds,segmentDelays,addZeros,scaleBy,resampFac,moveChunks,guess_best_scale1,isTimes)

% If like movie, scaleThisSignal=1
% else scaleThisSignal=0

signal=decimate(double(signal),decind);
if scaleThisSignal==1
    if isTimes==true
%         signal=resample(signal,floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);
%         temp=resample(signal,floor(scaleBy*resampFac),floor(resampFac));
        % different filter needed to prevent ringing artifacts and aliasing
        % Like movie
        normFc=0.98/max(floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);
        order=256*max(floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);
        beta=12; % smoothing
        lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
        lpFilt = lpFilt .* kaiser(order+1,beta)';
        lpFilt = lpFilt / sum(lpFilt);
        p=floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100; lpFilt = p * lpFilt;
        signal=resample(signal,floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100,lpFilt);
        % From initial alignment
        % signal=resample(signal,floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);
        % From final alignment
        %temp=resample(signal,floor(scaleBy*resampFac),floor(resampFac));
        normFc=0.98/max(floor(scaleBy*resampFac),floor(resampFac));
        order=256*max(floor(scaleBy*resampFac),floor(resampFac));
        beta=100; % smoothing
        lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
        lpFilt = lpFilt .* kaiser(order+1,beta)';
        lpFilt = lpFilt / sum(lpFilt);
        p=floor(scaleBy*resampFac); lpFilt = p * lpFilt;
        temp=resample(signal,floor(scaleBy*resampFac),floor(resampFac),lpFilt);
    else
        normFc=0.98/max(floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);
        order=256*max(floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);
        beta=12; % smoothing
        lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
        lpFilt = lpFilt .* kaiser(order+1,beta)';
        lpFilt = lpFilt / sum(lpFilt);
        p=floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100; lpFilt = p * lpFilt;

        signal=resample(signal,floor(mod(guess_best_scale1,1)*100)+floor((1*100)/100)*100,100);

        normFc=0.98/max(floor(scaleBy*resampFac),floor(resampFac));
        order=256*max(floor(scaleBy*resampFac),floor(resampFac));
        beta=5; % smoothing
        lpFilt = firls(order, [0 normFc normFc 1],[1 1 0 0]);
        lpFilt = lpFilt .* kaiser(order+1,beta)';
        lpFilt = lpFilt / sum(lpFilt);
        p=floor(scaleBy*resampFac); lpFilt = p * lpFilt;

        temp=resample(signal,floor(scaleBy*resampFac),floor(resampFac));
    end
    % cut off ringing artifact
    temp(end-10+1:end)=nan;
    signal=[nan(1,frontShift) temp];
    if movieToLength>length(signal)
        signal=[signal nan(1,movieToLength-length(signal))];
    end
else
    % Like arduino
    signal=[signal nan(1,shouldBeLength-length(signal))];
end

% [Xa,Ya] = alignsignals(X,Y)
% X is movie, Y is arduino
% If Y is delayed with respect to X, then D is positive and X is delayed by D samples.
% If Y is advanced with respect to X, then D is negative and Y is delayed by –D samples.

% disp('in second alignment');
% disp(length(signal))

outsignal=[];
firstInd=segmentInds(1);
outsignal=[outsignal zeros(1,firstInd-1)];
for i=1:length(segmentInds)-1
    currInd=segmentInds(i);
    if currInd+alignSegments-1>length(signal)
        currChunk=signal(currInd:end);
    else
        currChunk=signal(currInd:currInd+alignSegments-1);
    end
    if scaleThisSignal==1
        % Like movie
        if segmentDelays(i)>0
            % Delay is positive, so movie was shifted
            currChunk=[ones(1,segmentDelays(i))*currChunk(1) currChunk];
        else
            % Delay is negative, so arduino was shifted
        end
    else
        % Like arduino
        if segmentDelays(i)>0
            % Delay is positive, so movie was shifted
        else
            % Delay is negative, so arduino was shifted
            currChunk=[ones(1,-segmentDelays(i))*currChunk(1) currChunk];
        end 
    end
    currChunk=[currChunk currChunk(end)*ones(1,addZeros(i))];
    if moveChunks(i,2)>length(currChunk)
        currChunk=currChunk(moveChunks(i,1):end);
    else
        currChunk=currChunk(moveChunks(i,1):moveChunks(i,2));
    end
    outsignal=[outsignal currChunk];
end

end
