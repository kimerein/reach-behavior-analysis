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
endofDir=regexp(videoFile,'\');

%% Save discardFirstNFrames for rest of analysis
settings=autoReachAnalysisSettings(discardFirstNFrames);
save([videoFile(1:endofVfname(end)-1) '_autoReachSettings.mat'],'settings');

%% Get user-defined movie zones
setup_reach_coding(videoFile,discardFirstNFrames);
settings=setup_reach_coding_settings();
save([videoFile(1:endofVfname(end)-1) '_setupReachSettings.mat'],'settings');

%% Read data from user-defined zones
[out,zoneVals,reaches,pellets,eat,paw,fidget,settings]=extractEventsFromMovie([videoFile(1:endofVfname(end)-1) '_zones.mat'],videoFile,[]);

%% Re-format movie data as events
savehandles=reformatAutoClassifyOutput(out,zoneVals,reaches,pellets,eat,paw,fidget,settings);
save([videoFile(1:endofVfname(end)-1) '_savehandles.mat'],'savehandles');

%% Get Arduino output data
out=parseSerialOut_wrapper([videoFile(1:endofDir(end)) 'OUTPUT.txt'],[videoFile(1:endofVfname(end)-1) '_parsedOutput.mat']);
settings=arduinoSettings();
save([videoFile(1:endofVfname(end)-1) '_arduinoSettings.mat'],'settings');

%% Discard end of video
savehandles=discardLastNFrames(savehandles);

%% Align Arduino output data and data from video file
aligned=getAlignment(out,30,savehandles,[]);
save([videoFile(1:endofVfname(end)-1) '_aligned.mat'],'aligned');
settings=alignmentSettings();
save([videoFile(1:endofVfname(end)-1) '_alignmentSettings.mat'],'settings');
aligned2=aligned;
pause;

%% Clean up cue from movie
% aligned=cleanUpCue(aligned);
[aligned,cleanup]=cleanUpCue_basedOnArduino(aligned2);
save([videoFile(1:endofVfname(end)-1) '_cleanup_settings.mat'], 'cleanup');
pause;

%% Put Arduino output data and data from video file together
[status]=mkdir([videoFile(1:endofVfname(end)-1) '_processed_data']);
finaldata=integrateSDoutWithReaches(savehandles,out,30,aligned,[videoFile(1:endofVfname(end)-1) '_processed_data']);

%% Check for chewing of pellet (this should take a while -- pellet is large)
finaldata=checkForChewedPellet(finaldata);
alignment=finaldata;
save([videoFile(1:endofVfname(end)-1) '_processed_data' '/final_aligned_data.mat'],'alignment');

%% Plot results
tbt=plotCueTriggeredBehavior(finaldata,'cueZone_onVoff',0);
save([videoFile(1:endofVfname(end)-1) '_processed_data/tbt.mat'],'tbt');
settings=plotCueTriggered_settings();
save([videoFile(1:endofVfname(end)-1) '_plottingSettings.mat'],'settings');

end
