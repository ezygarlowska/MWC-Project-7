function [Sk,As] = Skewness_asymmetry(timeseries)
% This function calculated the skewness and asymmetry of for a
% given time-series of free-surface elevation (eta).  

% Detrend timeseries
eta = detrend(timeseries);

% Skewness
Sk = mean(eta.^3)./(mean(eta.^2).^(3/2));

% Asymmetry
As = mean(imag(hilbert(eta)).^3)./(mean(eta.^2).^(3/2));
end
