function analyzeReachVideo(videoFile,discardFirstNFrames)
% Note that matching Arduino output data file should be named 'OUTPUT.txt'
% and in same directory as videoFile

% Check the following files before running this code.
% These files determine the settings to be used for the analysis.
% 
% setup_reach_coding_settings.m  : get user-defined regions of video
% autoReachAnalysisSettings.m    : defining different types of events in movie
% alignmentSettings.m            : aligning movie events to Arduino output data
% arduinoSettings.m              : parse Arduino output data
% plotCueTriggered_settings.m    : how to plot experimental results

endofVfname=regexp(videoFile,'\.');
endofDir=regexp(videoFile,'/');

%% Save discardFirstNFrames for rest of analysis
autoReachAnalysisSettings(discardFirstNFrames);

%% Get user-defined movie zones
setup_reach_coding(videoFile,discardFirstNFrames);

%% Read data from user-defined zones
[out,zoneVals,reaches,pellets,eat,paw,fidget,settings]=extractEventsFromMovie([videoFile(1:endofVfname(end)-1) '_zones.mat'],videoFile,[]);

%% Re-format movie data as events
savehandles=reformatAutoClassifyOutput(out,zoneVals,reaches,pellets,eat,paw,fidget,settings);
save([videoFile(1:endofVfname(end)-1) 'savehandles.mat'],'savehandles');

%% Get Arduino output data
out=parseSerialOut_wrapper([videoFile(1:endofDir(end)) 'OUTPUT.txt'],[videoFile(1:endofVfname(end)) 'parsedOutput.mat']);

%% Align Arduino output data and data from video file
aligned=getAlignment(out,30,savehandles);

%% Clean up cue from movie
aligned=cleanUpCue(aligned);

%% Put Arduino output data and data from video file together    
[status]=mkdir([videoFile(1:endofVfname(end)-1) '_processed_data']);
finaldata=integrateSDoutWithReaches(savehandles,out,30,aligned,[videoFile(1:endofVfname(end)-1) '_processed_data']);  

%% Plot results
tbt=plotCueTriggeredBehavior(finaldata,'cue',1);
save([videoFile(1:endofVfname(end)-1) '_processed_data/tbt.mat'],'tbt');

end