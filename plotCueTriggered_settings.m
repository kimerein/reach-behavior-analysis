function settings=plotCueTriggered_settings()

% Inter-trial interval (ITI) settings for this experiment
% Will help with choosing one cue timepoint per trial
% in case of issues with aliasing of instantaneous cue
settings.maxITI=50; % in seconds, maximal ITI
settings.minITI=2; % in seconds, minimal ITI

% Min height of cue or pellet presented peak relative to neighboring
% samples
% Increase for more stringent cue / pellet presented detection
% Decrease to pick up more peaks
% settings.relativePeakHeight=1*10^35; 
settings.nStdDevs=1; % relative peak height in terms of standard deviations away from mean
% i.e., greater than nStdDevs * nanstd(data) away from mean

% Experiment block types:

% Get to choose one or more types of trial shading for the event plot
% shading_type can include
% 'ITI'     blocks determined by different ITI lengths
% settings.shading_type={'ITI'};
settings.shading_type={};

% For experiments with blocks of different trial lengths
% If 'ITI' is included in shading_type, will shade each block of a different
% trial length with a different  color
settings.blockITIThresh=[10]; % boundaries between ITI block lengths, in seconds
% e.g., if there are 2 trial lengths, blockITIThresh is the threshold that
% distinguishes short and long trials

% Shading colors
% The indices here correspond to indices into shading_type
settings.shading_colors={{[0.9 0.9 0.9],'none'}};                        

% How many time points to include from previous trial at the beginning of
% each trial
settings.pointsFromPreviousTrial=100;


% Plotting trial-by-trial averages:

% Which fields from data to plot as averaged trial-by-trial
% List the field names to plot (these are field names from data variable in
% plotCueTriggeredBehavior.m)
settings.plotfields={'cue', ...
                     'arduino_distractor',...
                     'pelletLoaded',...
                     'pelletPresented',...
                     'optoOn',...
                     'interlockSamples',...
                     'reachStarts',...
                     'reach_ongoing',...
                     'success_reachStarts',...
                     'drop_reachStarts',...
                     'miss_reachStarts',...
                     'pelletmissingreach_reachStarts',...
                     'eating'};

                 
% Plotting events:

% Which fields from data to plot as events
% Include field names in this list
% Will plot these events in this order, so event types later in the list
% will be plotted on top of the event types earlier in the list
settings.plotevents={'optoOn',...
                     'success_reachStarts_pawOnWheel',...
                     'drop_reachStarts_pawOnWheel',...
                     'miss_reachStarts_pawOnWheel',...
                     'pelletmissingreach_reachStarts',...
                     'pelletPresented',...
                     'success_reachStarts',...
                     'drop_reachStarts',...
                     'miss_reachStarts',...
                     'cue'};
                 
% Color for each event type
pawwheel_color='y'; % will fill in circle with this color if paw started from wheel for this reach
% Indices into eventColors and eventOutlines correspond to indices into plotevents
% eventColors gives fill of circle for each event type
% eventOutlines gives outline of circle for each event type
settings.eventColors={  'none',...
                        pawwheel_color,...
                        pawwheel_color,...
                        pawwheel_color...
                        [0.8 0.8 0.8],...
                        'k',...
                        'g',...
                        'r',...
                        'c',...
                        'b'};
settings.eventOutlines={'m',...
                        'g',...
                        'r',...
                        'c',...
                        [0.8 0.8 0.8],...
                        'none',...
                        'g',...
                        'none',...
                        'none',...
                        'b'};
                    
% Event thresholds for each event type
% All indices with data values above the threshold will be counted as events
% May want to specify different thresholds for each event type
% Indices into eventThresh correspond to indices into plotevents
settings.eventThresh={0.2,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.5,...
                      0.2,...
                      0.2,...
                      0.2,...
                      0.5};
                  
% Will only plot firstN events of each type in each trial
% For example, will only plot first moment of cue on in each trial
% Indices into firstN correspond to indices into plotevents
% If value is 'all', will plot all events of this type in trial
settings.firstN={'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 'all',...
                 1};
             
% Width of outline of event maker
settings.eventLineWidth=1.5;

% Some events, like the cue, are detected as the peak of the resampled cue
% data. The onset of the cue may precede this peak by a fixed time. All
% events will be shifted backwards in time by shiftBack seconds, in the figure. 
% For example, if shiftBack is 0.2 for plotevents 'cue', the cue event will be
% plotted 200 ms earlier. Indices into shiftBack correspond to indices 
% into plotevents. If in doubt, set shiftBack to zero.
settings.shiftBack={0,...
                    0,...
                    0,...
                    0,...
                    0,...
                    0.2,...
                    0,...
                    0,...
                    0,...
                    0.2};

                
% Plot histogram of events (or trial-by-trial average), with data types
% overlapping:
                
% Which fields from data to plot as overlapping averaged trial-by-trial 
% List the field names to plot (these are field names from data variable in
% plotCueTriggeredBehavior.m) in histoplotfields
% Note that will also shift each event time for each data type by 
% histoshiftBack
settings.histoplotfields={'cue',...
                          'pelletPresented',...
                          'reachStarts'};
settings.histoshiftBack={0.2,...
                         0.2,...
                         0};

end