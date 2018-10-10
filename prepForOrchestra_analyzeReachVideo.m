function prepForOrchestra_analyzeReachVideo(videoFile,discardFirstNFrames)

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