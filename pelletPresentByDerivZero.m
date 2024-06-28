function pellets=pelletPresentByDerivZero(pellets,zeroDerivRange,isOrchestra)

derivs=diff(pellets.rawData);
vals=pellets.rawData(1:end-1);
data=vals(derivs>zeroDerivRange(1) & derivs<zeroDerivRange(2));
data=data';

if ~isOrchestra
    % Plot the histogram
    figure;
    histogram(data, 50);
    title('Histogram of Data');
    xlabel('Data');
    ylabel('Frequency');
end

% Fit a Gaussian Mixture Model (GMM) to the data
gm = fitgmdist(data, 2);

% Extract the means and standard deviations of the two components
mu1 = gm.mu(1);
mu2 = gm.mu(2);
sigma1 = sqrt(gm.Sigma(1));
sigma2 = sqrt(gm.Sigma(2));

% Define a function for the Gaussian PDF
gaussian = @(x, mu, sigma) (1/(sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma).^2);

% Generate x values for plotting the Gaussian curves
x = linspace(min(data), max(data), 1000);

% Calculate the Gaussian curves
y1 = gaussian(x, mu1, sigma1) * gm.ComponentProportion(1);
y2 = gaussian(x, mu2, sigma2) * gm.ComponentProportion(2);

if ~isOrchestra
    % Plot the Gaussian curves on top of the histogram
    hold on;
    plot(x, y1 * length(data) * (x(2) - x(1)), 'r', 'LineWidth', 2);
    plot(x, y2 * length(data) * (x(2) - x(1)), 'b', 'LineWidth', 2);
end

% Find the threshold that best separates the two modes
threshold = fzero(@(x) gaussian(x, mu1, sigma1) - gaussian(x, mu2, sigma2), mean([mu1 mu2]));

% Plot the threshold on the histogram
%     yLimits = ylim;
%     plot([threshold threshold], yLimits, 'k--', 'LineWidth', 2);
%     hold off;

% Display the threshold
%     fprintf('Threshold: %.2f\n', threshold);

% Get the mean of the higher mode and use as cutoff for pellet present
% Note that this assumes the mean of the higher mode will be lower than the
% actual pellet present value, because there is some noise in between that
% the higher mode tries to capture
if mu1>mu2
    cutoff=(mu1+threshold)/2;
else
    cutoff=(mu2+threshold)/2;
end
if ~isOrchestra
    % Plot the cutoff on the histogram
    yLimits = ylim;
    plot([cutoff cutoff], yLimits, 'k--', 'LineWidth', 2);
end

% Use this cutoff to define pellet present or absent
pellets.gm=gm;
pellets.gm_cutoff=cutoff;
pellets.pelletPresent=pellets.rawData>cutoff;

pellets.fig=figure();
plot(pellets.rawData,'Color','r');
hold on;
plot(pellets.pelletPresent.*max(pellets.rawData,[],'all','omitnan'),'Color','g');

end