function fixReachVideo

% fixReachVideo.m

% this is a script

%% Settings

clear variables

videoFile='Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\20240820\X4\O2 output\2015-03-19 20-32-19-C.avi';
chronuxPath='C:\Users\sabatini\Documents\GitHub\chronux_2_11'; % path to Chronux
parsedOutputFile='Z:\MICROSCOPE\Kim\KER Behavior\By date\Low speed\20240820\X4\O2 output\2015-03-19 20-32-19-C_parsedOutput.mat';

%% Set up

% Set up video-specific variables
if ispc
    endofVfname=regexp(videoFile,'AVI');
    if isempty(endofVfname)
        endofVfname=regexp(videoFile,'avi');
    end
    endofVfname=endofVfname-1;
else
    endofVfname=regexp(videoFile,'AVI');
    if isempty(endofVfname)
        endofVfname=regexp(videoFile,'avi');
    end
    endofVfname=endofVfname-1;
end
endofDir=regexp(videoFile,sep);

%% Check O2 output STEP 1
% First, check reachesFig.fig
% Second, check pelletsFig.fig
% If either figure is problematic, run this section

% Variables to adjust:
discardMoreFramesAtBeginning=0; % Thow out this many more frames at the beginning of the video
chewThresh=1; % default is 1
vars.subtractExternalCue=true; % if external cue is lighting up whole field of view
vars.cueWasRamp=true; % if external cue was ramp, need to adjust cueZone_onVoff timing to catch BEGINNING rather than PEAK of ramp

if chewThresh~=1
    qans=questdlg('Chewing threshold is usually 1. Are you sure you want to proceed with a different value for chewThresh?');
    switch qans
        case 'Yes'
        case 'No'
            chewThresh=1;
        case 'Cancel'
            chewThresh=1;
    end
end

a=load([videoFile(1:endofVfname(end)-1) '_autoReachSettings.mat']);
settings=a.settings;
settings=autoReachAnalysisSettings(settings.discardFirstNFrames+discardMoreFramesAtBeginning,true,chronuxPath,chewThresh);
a=load([videoFile(1:endofVfname(end)-1) '_zoneVals.mat']);
zoneVals=a.zoneVals;
vars.videoFile=videoFile;
vars.endofVfname=endofVfname;
vars.endofDir=endofDir;
vars.zoneVals=zoneVals;
analyzeReachVideo_wrapper('extractEventsFromMovie',vars);

%% STEP 1.5
% Check whether reach accidentally grabbed distractor
% Check whether reach accidentally grabbed cue
% Check whether distractor accidentally grabbed cue (problem for external
% vis cue)

a=load([videoFile(1:endofVfname(end)-1) '_reaches.mat']);
reaches=a.reaches;
a=load([videoFile(1:endofVfname(end)-1) '_savehandles.mat']);
savehandles=a.savehandles;

figure();
plot(reaches.rawData-nanmin(reaches.rawData),'Color','k');
hold on; plot(((reaches.isReach-nanmin(reaches.isReach))./nanmax(reaches.isReach-nanmin(reaches.isReach))).*nanmax(reaches.rawData-nanmin(reaches.rawData)),'Color','r');
plot(savehandles.LEDvals-nanmin(savehandles.LEDvals),'Color','b');
legend({'black: raw reach data','red: classified as reach','blue: distractor'});

figure();
plot(reaches.rawData-nanmin(reaches.rawData),'Color','k');
hold on; plot(((reaches.isReach-nanmin(reaches.isReach))./nanmax(reaches.isReach-nanmin(reaches.isReach))).*nanmax(reaches.rawData-nanmin(reaches.rawData)),'Color','r');
plot(savehandles.cueZone-nanmin(savehandles.cueZone),'Color','b');
legend({'black: raw reach data','red: classified as reach','blue: cue'});

figure();
plot(savehandles.LEDvals-nanmin(savehandles.LEDvals),'Color','b');
hold on;
plot(savehandles.cueZone-nanmin(savehandles.cueZone),'Color','r');
legend({'blue: distractor','red: cue'});

figure(); plot(reaches.rawData,'Color','r'); hold on; plot(savehandles.LEDvals,'Color','b');
hold on; plot(savehandles.pelletPresent,'Color','g');
hold on; plot(savehandles.pelletPresent*8e5,'Color','g');
hold on; plot(savehandles.optoZone,'Color','k');
hold on; plot(savehandles.cueZone,'Color','m');
%% STEP 2 -- Alignment of Arduino and Movie
% If there are problems with the alignment (or if you've re-run STEP 1), run this

% Variables to adjust:
discardFramesAtEnd=0; % Thow out this many frames at the end of the video

% Threshold for distinguishing LED distractor on vs off
% The threshold will be min(LED distractor) + fractionRange*range(LED
% distractor)
fractionRange=0.5;

% The following values help the alignment by giving an estimate of when
% the movie fits into the arduino data.
% If isInSecondHalf is false, the code will start looking for an appropriate
% alignment of the movie data at the beginning of the arduino data.
% If, in fact, the movie comes in the second half of the arduino data
% stream, indicate this by setting isInSecondHalf to true.
isInSecondHalf=false; % set this to true if movie matches a later section of arduino data stream

% For fractionThroughArduino ...
% Where in the arduino data stream does the movie begin?
% If isInSecondHalf is 1, we will discard the first nth of the arduino data
% stream before looking for an alignment with the movie, where n is
% fractionThroughArduino.
% This helps code find the correct alignment.
% For example, if the movie begins 75% of the way through the arduino data
% stream, set fractionThroughArduino to 3/4.
fractionThroughArduino=0.6;

% The code will try different scalings of the movie data onto the arduino
% data. An initial guess at the correct scaling will be chosen based on a
% preliminary alignment. The code will then try to further refine this
% estimate (called guess_best_scale) by trying different scalings similar
% to this best guess. The code will try all scalings between
% tryscales=guess_best_scale+try_scale1:tryinc:guess_best_scale+try_scale2
tryinc=0.00005; % this is the increment for trying different scalings of movie onto arduino data
try_scale1=0; % minimum scaling to try (+guess_best_scale)
try_scale2=0.02; % maximum scaling to try (+guess_best_scale)

% Change delays to try
try_delay1=0;
try_delay2=300;


if ~isempty(parsedOutputFile)
    a=load(parsedOutputFile);
    out=a.out;
else
    a=load([videoFile(1:endofVfname(end)-1) '_parsedOutput.mat']);
    out=a.out;
end
a=load([videoFile(1:endofVfname(end)-1) '_savehandles.mat']);
savehandles=a.savehandles;
if discardFramesAtEnd>0
    savehandles.discardLastN=discardFramesAtEnd;
end

vars.savehandles=savehandles;
vars.out=out;
vars.videoFile=videoFile;
vars.endofVfname=endofVfname;
vars.fractionRange=fractionRange;
vars.isInSecondHalf=isInSecondHalf;
vars.fractionThroughArduino=fractionThroughArduino;
vars.tryinc=tryinc;
vars.try_scale1=try_scale1;
vars.try_scale2=try_scale2;
vars.try_delay1=try_delay1;
vars.try_delay2=try_delay2;

analyzeReachVideo_wrapper('getAlignment',vars);

%% STEP 3 -- Find cue events in movie
% If there are problems with cue detection (or if you've re-run STEP 1 or
% STEP 2), run this

% Variables to adjust:
minProm=20000; % minimum height of "cue" signal in movie

a=load([videoFile(1:endofVfname(end)-1) '_aligned.mat']);
aligned2=a.aligned;

vars.minProm=minProm;
vars.aligned2=aligned2;
vars.videoFile=videoFile;
vars.endofVfname=endofVfname;

if vars.subtractExternalCue==true
    % cue and distractor never on at the same time
    aligned2.backup_cueZone=aligned2.cueZone;
    aligned2.temp_movie_distractor=aligned2.movie_distractor;
    figure(); plot(aligned2.cueZone,'Color','b');
    base=input('Cue baseline: ');
    % Expand movie_distractor by a few inds
    resampleSlop=1000; % in ms
    resampleSlopInds=ceil((resampleSlop/1000)/0.03);
    f=aligned2.temp_movie_distractor>0.5;
    for i=1+resampleSlopInds:length(f)-1-resampleSlopInds
        if f(i)==0 & f(i+1)==1
            aligned2.temp_movie_distractor(i-resampleSlopInds:i+resampleSlopInds)=1;
        elseif f(i)==1 & f(i+1)==0
            aligned2.temp_movie_distractor(i-resampleSlopInds:i+resampleSlopInds)=1;
        end
    end
    aligned2.cueZone(aligned2.temp_movie_distractor>0.5)=base;
    vars.aligned2=aligned2;
end
    
analyzeReachVideo_wrapper('getCue',vars);

if vars.subtractExternalCue==true && vars.cueWasRamp==true
    analyzeReachVideo_wrapper('fixRampExternalCue',vars);
end

%% STEP 4 -- Make final organized data structure in folder

if ~isempty(parsedOutputFile)
    a=load(parsedOutputFile);
    out=a.out;
else
    a=load([videoFile(1:endofVfname(end)-1) '_parsedOutput.mat']);
    out=a.out;
end
a=load([videoFile(1:endofVfname(end)-1) '_savehandles.mat']);
savehandles=a.savehandles;
a=load([videoFile(1:endofVfname(end)-1) '_aligned.mat']);
aligned=a.aligned;

vars.out=out;
vars.savehandles=savehandles;
vars.aligned=aligned;

analyzeReachVideo_wrapper('organizeData',vars);

%% STEP 5 -- Run SVM to do final classification of success versus drop

doThisManually=true;

addpath(genpath(chronuxPath));
if doThisManually==true
    fixDropVSuccess([videoFile(1:endofVfname(end)-1) '_processed_data'],[videoFile(1:endofVfname(end)-1)],[],false);
else
    fixDropVSuccess([videoFile(1:endofVfname(end)-1) '_processed_data'],[videoFile(1:endofVfname(end)-1)],[]);
end
fid=fopen([videoFile(1:endofVfname(end)-1) '_processed_data\fixed_miss_v_grab.txt'],'wt');
fclose(fid);

%% RESET PATH

rmpath(genpath(chronuxPath));

