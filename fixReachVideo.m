% fixReachVideo.m

% this is a script

%% Settings

clear variables

videoFile='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior\By date\Low speed\20200218\Nov_ON\O2 output\2020-02-18 15-51-49-C.AVI';
chronuxPath='C:/Users/kim/Documents/MATLAB/chronux_2_11'; % path to Chronux
parsedOutputFile='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\KER Behavior\By date\Low speed\20200218\Nov_ON\O2 output\2020-02-18 15-51-49-C_parsedOutput.mat';

%% Set up 

% Set up video-specific variables
endofVfname=regexp(videoFile,'\.');
endofDir=regexp(videoFile,'\');

%% Check O2 output STEP 1
% First, check reachesFig.fig
% Second, check pelletsFig.fig
% If either figure is problematic, run this section

% Variables to adjust:
discardMoreFramesAtBeginning=0; % Thow out this many more frames at the beginning of the video

a=load([videoFile(1:endofVfname(end)-1) '_autoReachSettings.mat']);
settings=a.settings;
settings=autoReachAnalysisSettings(settings.discardFirstNFrames+discardMoreFramesAtBeginning,true,chronuxPath);
a=load([videoFile(1:endofVfname(end)-1) '_zoneVals.mat']);
zoneVals=a.zoneVals;
vars.videoFile=videoFile;
vars.endofVfname=endofVfname;
vars.endofDir=endofDir;
vars.zoneVals=zoneVals;
analyzeReachVideo_wrapper('extractEventsFromMovie',vars);

%% STEP 2 -- Alignment of Arduino and Movie
% If there are problems with the alignment (or if you've re-run STEP 1), run this

% Variables to adjust:
discardFramesAtEnd=0; % Thow out this many frames at the end of the video

% Threshold for distinguishing LED distractor on vs off
% The threshold will be min(LED distractor) + fractionRange*range(LED
% distractor)
fractionRange=0.3;

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

analyzeReachVideo_wrapper('getAlignment',vars);

%% STEP 3 -- Find cue events in movie
% If there are problems with cue detection (or if you've re-run STEP 1 or
% STEP 2), run this

% Variables to adjust:
minProm=100000; % minimum height of "cue" signal in movie

a=load([videoFile(1:endofVfname(end)-1) '_aligned.mat']);
aligned2=a.aligned;

vars.minProm=minProm;
vars.aligned2=aligned2;
vars.videoFile=videoFile;
vars.endofVfname=endofVfname;

analyzeReachVideo_wrapper('getCue',vars);

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
