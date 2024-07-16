function [chi2stat] = AAA_test_unif(data, fs)

% Number of bins
numBins = round(max(data)/fs);

if numBins == 0
    numBins = round(length(data)/2);
end

% Calculate observed frequencies
[observedFrequencies, edges] = histcounts(data, numBins);

% Calculate expected frequencies
expectedFrequencies = numel(data) / numBins * ones(1, numBins);

% Compute the Chi-Square statistic
chi2stat = sum((observedFrequencies - expectedFrequencies).^2 ./ expectedFrequencies);

% Degrees of freedom
df = numBins - 1;

% Critical value from Chi-Square distribution
alpha = 0.05;  % Significance level
criticalValue = chi2inv(1 - alpha, df);

% Display results
% fprintf('Chi-Square Statistic: %.4f\n', chi2stat);
% fprintf('Critical Value: %.4f\n', criticalValue);
% 
% if chi2stat > criticalValue
%     disp('Reject the null hypothesis: The data does not follow a uniform distribution.');
% else
%     disp('Fail to reject the null hypothesis: The data follows a uniform distribution.');
% end
