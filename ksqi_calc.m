function ksqi=ksqi_calc(signal)
% Function to calculate Kurtosis Based SQI
% Input signal = signal in
% Output fsqi = Derived SQI

% Find mean of signal
mu = mean(signal);

% Alpha Parameter for IIR Filter
alpha = 0.05;

% Numerator of Kurtosis
num = filter(alpha,[1,alpha-1],(signal-mu).^4);

% Denominator of Kurtosis
den = filter(alpha,[1,alpha-1],(signal-mu).^2);

% Calculate KSQI 
ksqi = (num./den.^2);
ksqi(ksqi<0) = 0;

end

