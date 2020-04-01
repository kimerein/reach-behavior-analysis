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

%% Write these movie zones to a text file
% (will use OpenCV in Python to read in video file)
a=load([videoFile(1:endofVfname(end)-1) '_zones.mat']);
zones=a.zones;
fid=fopen([videoFile(1:endofVfname(end)-1) '_zonesAsText.txt'],'wt');
for i=1:length(zones)
    takeImageValue=zones(i).takeImageValue;
    isin=zones(i).isin;
    dim1points=zones(i).dim1points;
    dim2points=zones(i).dim2points;
    dim1coords=dim1points(isin);
    dim2coords=dim2points(isin);
    tuplesString='[';
    for j=1:length(dim1coords)
        thisTuple=['[' num2str(dim1coords(j),'%i') ',' num2str(dim2coords(j),'%i') ']'];
        if j==length(dim1coords)
            tuplesString=[tuplesString thisTuple ']'];
        else
            tuplesString=[tuplesString [thisTuple ',']];
        end
    end
    fprintf(fid,'%s %s\n\n',takeImageValue,tuplesString);
end
fclose(fid);   
