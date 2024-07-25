function aligned=fixRampExternalCue(aligned)

[modeDelay,peaks,delays,rampBeginnings] = calculateModeDelay(aligned.cueZone, find(aligned.cueZone_onVoff==1));
% Make rampBeginnings cueZone_onVoff
figure(); plot(aligned.cueZone_onVoff,'Color','b'); hold on;
aligned.cueZone_onVoff(:)=0;
aligned.cueZone_onVoff(rampBeginnings)=1;
plot(aligned.cueZone_onVoff,'Color','r');

end

function rampBeginnings = findRampBeginnings(values, peaks)
    % Preallocate the rampBeginnings array
    rampBeginnings = NaN(size(peaks));
    
    % Loop through each peak to find the beginning of the ramp
    for i = 1:length(peaks)
        peakIndex = peaks(i);
        
        % Initialize the ramp beginning index
        rampBeginIndex = peakIndex;
        
        % Trace back to find the start of the ramp
        for j = peakIndex-1:-1:1
            if values(j) < values(j+1)
                rampBeginIndex = j;
            else
                break; % Stop if a decrease is found
            end
        end
        
        % Store the ramp beginning index
        rampBeginnings(i) = rampBeginIndex;
    end
end

function [modeDelay,peaks,delays,rampBeginnings] = calculateModeDelay(values, indices)
    % Find the nearest peaks for the given indices
    peaks = findNearestPeaks(values, indices);
    
    % Find the ramp beginnings for each peak
    rampBeginnings = findRampBeginnings(values, peaks);
    
    % Calculate the delay between ramp beginnings and peaks
    delays = peaks - rampBeginnings;
    
    % Calculate the mode of the delays
    modeDelay = mode(delays);
end

function nearestPeaks = findNearestPeaks(values, indices)
    % Preallocate the nearestPeaks array with NaNs
    nearestPeaks = NaN(size(indices));
    
    % Loop through each index in the indices array
    for i = 1:length(indices)
        idx = indices(i);
        
        % Initialize the peak value as the value at the current index
        peakValue = values(idx);
        
        % Initialize the peak index
        peakIndex = idx;
        
        % Check for the nearest peak after the current index
        for j = idx:length(values)
            if values(j) >= peakValue
                peakValue = values(j);
                peakIndex = j;
            else
                break; % Stop if a decrease is found
            end
        end
        
        % Store the nearest peak index
        nearestPeaks(i) = peakIndex;
    end
end
