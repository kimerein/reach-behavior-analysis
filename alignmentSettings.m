function settings=alignmentSettings()

% Settings for getAlignment.m and discardLastNFrames.m
% Note that microSD (Arduino) output is timed in ms
% Whereas video is timed in frames per sec

% Discard the last N frames of the movie where N is discardLastN
settings.discardLastN=0;

% Threshold for distinguishing LED distractor on vs off
% The threshold will be min(LED distractor) + fractionRange*range(LED
% distractor)
settings.fractionRange=0.75;

% Minimum time between distractor LED on intervals
settings.minLEDinterval=1; % in seconds
settings.minCueInterval=0.2; % in seconds

% Need to resample movie and arduino data so that indices represent
% matching times
settings.arduino_fs=1000; % arduino data sampling rate in Hz
settings.movie_fs=30; % movie data sampling rate in Hz
settings.scale_factor=floor(settings.arduino_fs/settings.movie_fs);

% Choose the following two numbers based on approximate relationship between
% sampling rate of arduino data and sampling rate of movie data.
% For example, if movie rate is 30 frames per second and arduino data is
% timed in ms, 1000 ms/30 ms = 33 ... scale factor is 33.
settings.arduino_dec=settings.scale_factor;
settings.movie_dec=1;

% Throw out distractor LED durations in movie or arduino less than this
% many ms
<<<<<<< Updated upstream
% settings.useDistractorThresh=167; % in ms
settings.useDistractorThresh=175; % in ms
% settings.useDistractorThresh=160; % in ms
=======
settings.useDistractorThresh=200; % in ms
% settings.useDistractorThresh=175; % in ms
% settings.useDistractorThresh=165; % in ms
>>>>>>> Stashed changes
% settings.useDistractorThresh=150; % in ms

% If, for example, experimenter forget to include LED in movie frame at beginning of
% experiment, but Arduino was on, need to discard beginning of arduino
% distractor LED. Discard this much time from the beginning of arduino
% distractor LED.
% settings.discardTimeArduinoLED=3; % in seconds
settings.discardTimeArduinoLED=0; % in seconds

% Maxlag for initial alignment
% May help constrain alignment for better results
% This is the maximum window size to find the estimated delay, D, between
% the two input signals, arduino_LED and movie_LED
% maxlag from alignsignals Matlab function
% Set this to empty array, [], if don't want to constrain initial alignment
settings.maxlagForInitialAlign=50;

% The following values help the alignment by giving an estimate of when
% the movie fits into the arduino data.
% If isInSecondHalf is false, the code will start looking for an appropriate
% alignment of the movie data at the beginning of the arduino data.
% If, in fact, the movie comes in the second half of the arduino data
% stream, indicate this by setting isInSecondHalf to true.
settings.isInSecondHalf=false; % set this to true if movie matches a later section of arduino data stream

% For fractionThroughArduino ...
% Where in the arduino data stream does the movie begin? 
% If isInSecondHalf is 1, we will discard the first nth of the arduino data
% stream before looking for an alignment with the movie, where n is
% fractionThroughArduino.
% This helps code find the correct alignment.
% For example, if the movie begins 75% of the way through the arduino data
% stream, set fractionThroughArduino to 3/4.
settings.fractionThroughArduino=0.1; 

% The code will try different scalings of the movie data onto the arduino
% data. An initial guess at the correct scaling will be chosen based on a
% preliminary alignment. The code will then try to further refine this
% estimate (called guess_best_scale) by trying different scalings similar 
% to this best guess. The code will try all scalings between
% tryscales=guess_best_scale+try_scale1:tryinc:guess_best_scale+try_scale2
settings.tryinc=0.00005; % this is the increment for trying different scalings of movie onto arduino data
settings.try_scale1=-0.01;
settings.try_scale2=0.01;  
% If the preliminary alignment seems to produce an under-scaling of movie
% data with respect to arduino data, increase try_scale1 and try_scale2.
% If the preliminary alignment seems to produce an over-scaling of movie
% data with respect to arduino data, decrease try_scale1 and try_scale2.

% Similarly, the code will try different delays of the movie data with
% respect to the arduino data. An initial guess at the correct delay is
% chosen based on a preliminary alignment. The code will then try to
% further refine this estimate (called guess_best_delay) by trying
% different delays similar to this best guess. The code will try all delays
% between
% trydelays=guess_best_delay+try_delay1:guess_best_delay+try_delay2;
settings.try_delay1=-150;
settings.try_delay2=150;

% The movie DVR occasionally skips. For final alignment, code will subtly 
% shift sub-sections of movie data to better match arduino data 
% settings.alignSegments=600; % how many indices in each sub-section to independently align
settings.alignSegments=1750; % how many indices in each sub-section to independently align
% settings.alignSegments=2500; % how many indices in each sub-section to independently align
% For more precise local alignment, decrease alignSegments. For more
% precise global alignment, increase alignSegments.

% The code will automatically align the distractor and cue
% Here we specify additional data to align
% These data types are from arduino
% Align all fields with fieldname alignField(i).name
% Any field named 'cue' or 'falseCueOn' will be automatically aligned (no
% need to include the cue in this list)
settings.alignField(1).name='pelletLoaded';
settings.alignField(1).fromarduino=1;
settings.alignField(2).name='pelletPresented';
settings.alignField(2).fromarduino=1;
settings.alignField(3).name='optoOn';
settings.alignField(3).fromarduino=1;
settings.alignField(4).name='interlockSamples';
settings.alignField(4).fromarduino=1;
settings.alignField(5).name='solenoidOn';
settings.alignField(5).fromarduino=1;
% These data types are from movie
% Note that the names for these should match the names of fields in
% zoneVals
settings.alignField(6).name='optoZone';
settings.alignField(6).fromarduino=0;
settings.alignField(7).name='lickZone';
settings.alignField(7).fromarduino=0;
settings.alignField(8).name='cueZone';
settings.alignField(8).fromarduino=0;
settings.alignField(9).name='isGrooming';
settings.alignField(9).fromarduino=0;
settings.alignField(10).name='isChewing';
settings.alignField(10).fromarduino=0;

% Whether to add back all LED distractor events (including short ones, on
% the edge of what movie frame rate can capture) before final alignment
settings.alignWithAllEvents=1;

% Whether to allow user to manually zero out LED distractor on intervals
settings.doManualZeroOut=0; % 1 if want this user option, 0 otherwise

% Whether to use LED distractor to further align movieframeinds
settings.donotdoalign=0; % 0 if want further alignment of movieframeinds based on LED distractor, else 1

% How many frames does it take for LED to turn on?
% This will depend on the movie frame rate
settings.maxNFramesForLEDtoChange=2; % LED takes, at most, this many frames to turn on in movie

end