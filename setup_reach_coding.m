function [zones,custom_answers]=setup_reach_coding(varargin)

% Gets user-defined zones of movie
% See setup_reach_coding_settings.m for which zones to get
%
% varargin{1} is movie file
% varargin{2} is the number of frames that will be discarded from beginning
% of movie; use this to choose a time window in the movie for which all
% zones are well-defined

% Define some global variables for communication with zoneGUI
global continueAnalysis
continueAnalysis=0;

% Settings
settings=setup_reach_coding_settings();

% Get file name of video with reaches
filename=varargin{1};
discardFirstNFrames=varargin{2};

if isempty(discardFirstNFrames)
    error('Second argument discardFirstNFrames should be an integer');
end

% Instructions to user
continuebutton=questdlg(settings.prompt1,'Instructions 1','Yes','Cancel','Cancel');
switch continuebutton
    case 'Yes'
    case 'Cancel'
        return
end

% If filename indicates just a snapshot, then load image instead of video
if ~isempty(regexp(filename,'.png'))
    currentFrameNumber=1;
    allframes(:,:,:,currentFrameNumber)=imread(filename);
    custom_answers=[];
else
    % Read beginning of movie and discard unused frames
    % videoFReader = vision.VideoFileReader(filename,'PlayCount',1,'ImageColorSpace','YCbCr 4:2:2');
    videoFReader = vision.VideoFileReader(filename,'PlayCount',1);
    n=discardFirstNFrames; % How many frames to read initially
    if discardFirstNFrames>0
        for i=1:n
            [~,EOF]=step(videoFReader);
            if EOF==true
                break
            end
        end
    end

    n=settings.framesPerChunk;
    for i=1:settings.framesPerChunk % How many frames to read initially
        [frame,EOF]=step(videoFReader);
        if EOF==true
            allframes=allframes(:,:,:,i-1);
            break
        end
        if i==1
            allframes=nan([size(frame,1) size(frame,2) size(frame,3) n]);
        end
        allframes(:,:,:,i)=frame;
    end

    % Play movie until at good frame for zone definitions
    fig=implay(allframes,30);
    fig.Parent.Position=[100 100 800 800];
    pause;
    if ~isempty(regexp(version,'2017b','once')) || ~isempty(regexp(version,'2018b','once')) ||  ~isempty(regexp(version,'2020b','once')) || ~isempty(regexp(version,'2021a','once')) || ~isempty(regexp(version,'2021b','once'))
        currentFrameNumber=fig.DataSource.Controls.CurrentFrame;
    else
        currentFrameNumber=fig.data.Controls.CurrentFrame;
    end

    % Instructions to user
    continuebutton=questdlg(settings.prompt2,'Instructions 1','Yes','Cancel','Cancel');
    switch continuebutton
        case 'Yes'
        case 'Cancel'
            return
    end

    % Ask user other questions
    custom_answers=[];
    if isfield(settings,'n_custom')
        if settings.n_custom>0
            custom_answers=nan(1,settings.n_custom);
            for i=1:settings.n_custom
                continuebutton=questdlg(settings.prompts_custom{i},['Question ' num2str(i)],'Yes','No','No');
                switch continuebutton
                    case 'Yes'
                        custom_answers(i)=1;
                    case 'No'
                        custom_answers(i)=0;
                end
            end
        end
    end

    % Close implay fig, reopen an image so user can draw in zones
    close(fig);
end

zones=settings.zones;
for i=1:length(zones)
    disp(['Getting ' zones(i).name]);
    currZone=getZone(zones(i).prompt, zones(i).name, allframes(:,:,:,currentFrameNumber));
    [cols,rows]=find(ones(size(allframes,1),size(allframes,2))>0);
    zones(i).isin=inpolygon(rows,cols,currZone(:,1),currZone(:,2));
    zones(i).dim1points=rows;
    zones(i).dim2points=cols;
end

% Save zones
if settings.save==1
    endoffname=regexp(filename,'\.');
    save([filename(1:endoffname(end)-1) '_zones.mat'],'zones');
end

function zone=getZone(prompt1, prompt2, frame)

global zoneVertices
global continueAnalysis

zoneFig=perchZoneGUI(frame,[prompt1 ' Press "Done" after have defined vertices.']);
disp(['Press "enter" once have defined ' prompt2 '.']);
pause;
zone=zoneVertices;

if continueAnalysis==1
    disp([prompt2 ' succesfully defined']);
else
    disp(['failed to define ' prompt2]);
end

close(zoneFig);