function fsqi = fsqi_calc(signal,Fs)
% Function to calculate Frequency Based SQI of a given signal using Power
% in frequnecy bands
% Input signal = signal in
% Input Fs = Sample Rate
% Output fsqi = Derived SQI

alpha = 0.05; % parameter for IIR Filter

% Bandpass Filters to get signal in bands
relevant = bandpass(signal,[5 15],Fs);
high = bandpass(signal,[5 50],Fs);

% Square to get power
power_r = relevant.^2;
power_h = high.^2;

% IIR filter to average Power for per Sample SQI
num = filter(alpha,[1,alpha-1],power_r);
den = filter(alpha,[1,alpha-1],(power_h));

% Output is ratio of band (5 - 15) to high (5 - 50)
fsqi = (num./den);

end