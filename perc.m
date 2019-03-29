function [upper,lower,rel] = perc(x_hat,x,N,percentile)
% Calculate percentiles of errors
% Input x_hat = estimate 
% Input x = actual signal
% Input N = number of intervals
% Input percentile = percentile to calculate
% Output upper = upper perccentile
% Output lower = lower percentile
% Output rel = relative error

delta = sort((x_hat-x)/N);
delta = delta(~isnan(delta));
percentile = round(length(delta)*((100-percentile)/100));
delta = delta(percentile:end-percentile);
upper = max(delta);
lower = min(delta);
rel = mean(delta);

end

