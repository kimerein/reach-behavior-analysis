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

maxFrames=50000; % a single video file should not contain more video frames than this 

addpath(genpath('reachBehavior'));
addpath(genpath('reach_behavior_analysis'));
disp('added to path');

% Load data
f=dir(pwd);
j=1;
movie_names=[];
for i=1:length(f)
    currname=f(i).name;
    if ~isempty(regexp(currname,'\.AVI', 'once'))
        % is movie file
        movie_names{j}=currname;
        j=j+1;
    end
end

% For each movie file ... process with default settings
for i=1:length(movie_names)
    videoFile=movie_names{i};
    videoFile=[pwd '/' videoFile];
    disp(videoFile);
    endofVfname=regexp(videoFile,'\.');
    endofDir=regexp(videoFile,'/');
    
    % Read in autoReachAnalysisSettings
    a=load([videoFile(1:endofVfname(end)-1) '_autoReachSettings.mat']);
    settings=autoReachAnalysisSettings(a.settings.discardFirstNFrames);
    
    % Load names of zones
    a=load([videoFile(1:endofVfname(end)-1) '_zones.mat']);
    zones=a.zones;

    % Load data from user-defined movie zones (note that video data now
    % extracted by OpenCV/Python, prior to running this script)
    % Open file
    fid=fopen([videoFile(1:endofVfname(end)-1) '_zonesForMatlab.txt']);
    cline=fgetl(fid);
    dataForZoneVals=nan(maxFrames,length(zones));
    j=1;
    while cline~=-1
        % is -1 at eof
        a=eval(['[' cline ']']);
        if isempty(a)
            break
        end
        dataForZoneVals(j,:)=a;
        j=j+1;
        cline=fgetl(fid);
    end
    for j=1:size(dataForZoneVals,2)
        temp=dataForZoneVals(:,j)';
        temp(temp==-10)=nan;
        zoneVals.(zones(j).analysisField)=temp;
    end
    save([videoFile(1:endofVfname(end)-1) '_zoneVals.mat'],'zoneVals');
        
    % Extract features from video zones
    [out,zoneVals,reaches,pellets,eat,paw,fidget,settings]=extractEventsFromMovie([videoFile(1:endofVfname(end)-1) '_zones.mat'],videoFile,zoneVals);
    
    % Re-format movie data as events
    savehandles=reformatAutoClassifyOutput(out,zoneVals,reaches,pellets,eat,paw,fidget,settings);
    save([videoFile(1:endofVfname(end)-1) '_savehandles.mat'],'savehandles');
    
    % Get Arduino output data
    out=parseSerialOut_wrapper([videoFile(1:endofVfname(end)-1) '_OUTPUT.TXT'],[videoFile(1:endofVfname(end)-1) '_parsedOutput.mat']);
    settings=arduinoSettings();
    save([videoFile(1:endofVfname(end)-1) '_arduinoSettings.mat'],'settings');
    
    % Discard end of video
    savehandles=discardLastNFrames(savehandles);
    
    % Align Arduino output data and data from video file
    aligned=getAlignment(out,30,savehandles,[]);
    save([videoFile(1:endofVfname(end)-1) '_aligned.mat'],'aligned');
    settings=alignmentSettings();
    save([videoFile(1:endofVfname(end)-1) '_alignmentSettings.mat'],'settings');
    aligned2=aligned;
    
    % Clean up cue from movie
    % aligned=cleanUpCue(aligned);
    [aligned,cleanup]=cleanUpCue_basedOnArduino(aligned2);
    save([videoFile(1:endofVfname(end)-1) '_cleanup_settings.mat'], 'cleanup');
    close all
    
    % Put Arduino output data and data from video file together
    [status]=mkdir([videoFile(1:endofVfname(end)-1) '_processed_data']);
    finaldata=integrateSDoutWithReaches(savehandles,out,30,aligned,[videoFile(1:endofVfname(end)-1) '_processed_data']);
    
    % Check for chewing of pellet (this should take a while -- pellet is large)
    finaldata=checkForChewedPellet(finaldata);
    alignment=finaldata;
    save([videoFile(1:endofVfname(end)-1) '_processed_data' '/final_aligned_data.mat'],'alignment');
    
    % Plot results
    if ~isfield(finaldata,'isGrooming')
        finaldata.isGrooming=zeros(size(finaldata.cue));
    end
    tbt=plotCueTriggeredBehavior(finaldata,'cueZone_onVoff',0);
    save([videoFile(1:endofVfname(end)-1) '_processed_data/tbt.mat'],'tbt');
    settings=plotCueTriggered_settings();
    save([videoFile(1:endofVfname(end)-1) '_plottingSettings.mat'],'settings');
    close all
    
end

