function result_signal=removeExternalCue(signal1, signal2, signalMask, lpfilter)

% Sample signals (replace with actual signals)
% signal1 = randn(1, 100); % Example signal 1
% signal2 = 2 * signal1 + randn(1, 100); % Example signal 2 (contains signal 1 plus noise)

% Get info from user
figure(); plot(signalMask); title('signal mask');
maskThresh=input('Thresh: ');
figure(); plot(signal2); title('signal 2');
signal2Baseline=input('Signal 2 baseline: ');
signal2touse=signal2.*(signalMask>maskThresh);
signal2touse(signalMask<maskThresh)=signal2Baseline;
signal2touse=signal2touse-signal2Baseline;
figure(); plot(signal2touse); title('signal 2 masked and baseline-subtracted');
figure(); plot(signal1); title('signal 1');
signal1Baseline=input('Signal 1 baseline: ');
tosubtract=signal1-signal1Baseline;

% Define the range of scaling constants to try
scaling_constants = linspace(0, 2, 10000); % Adjust range and number of points as needed

% Initialize variables to store the best scaling constant and the minimal residual
best_scaling_constant = 0;
minimal_residual = inf;

% Iterate through the range of scaling constants
for scaling_constant = scaling_constants
    % Scale signal 1
    scaled_signal1 = scaling_constant * tosubtract;
    
    % Calculate the residual
    residual = abs(sum(signal2touse - scaled_signal1,'all','omitnan'));
    
    % Negative penalty
    % Penalize negative-going signal
%     negpen=sum(signal2touse - scaled_signal1<0,'all','omitnan');
%     if negpen>5000
%         continue
%     end
    
    % Update the best scaling constant if the current residual is smaller
    if residual < minimal_residual
        minimal_residual = residual;
        best_scaling_constant = scaling_constant;
    end
end

% Scale signal 1 with the best scaling constant
scaled_signal1 = best_scaling_constant * tosubtract;

% Subtract scaled signal 1 from signal 2
result_signal = signal2 - scaled_signal1;

if lpfilter==true
    cutoff_freq=3;
    fs=33;
    filter_order=4;
    result_signal(isnan(result_signal))=mean(result_signal,'all','omitnan');
    [b,a]=butter(filter_order, cutoff_freq/(fs/2), 'low');
    result_signal2=filtfilt(b,a,result_signal);
    figure(); plot(result_signal,'Color','k'); hold on; plot(result_signal2,'Color','r');
    pelletThresh=input('Thresh for pellet: ');
    result_signal=result_signal2;
end

% Plot the signals for visualization
figure;
subplot(3, 1, 1);
plot(signal1);
title('Signal 1');
subplot(3, 1, 2);
plot(signal2);
title('Signal 2');
subplot(3, 1, 3);
plot(result_signal);
title('Result Signal (Signal 2 minus scaled Signal 1)');
hold on;
plot((result_signal>pelletThresh).*nanmax(result_signal),'Color','g');

% Display best scaling constant and minimal residual
disp(['Best scaling constant: ', num2str(best_scaling_constant)]);
disp(['Minimal residual: ', num2str(minimal_residual)]);

result_signal=result_signal>pelletThresh;
end