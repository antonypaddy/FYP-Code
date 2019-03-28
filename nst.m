function [out,tm] = nst(signal,noise,ann,Fs,SNR,record_name)
% Function to add Noise to Signal using guide from PhysioNet 
% Input signal = signal to add noise to 
% Input noise = noise file 
% Input ann = annotations
% Input Fs = Sample Rate
% Input SNR = level of noise in dB
% Input record_name = string containing an extension to the file name
% Output out = signal out
% Output tm = time vector

% Size we are taking samples of signal power from
window = round(Fs*0.1);

out = zeros(size(signal));
new = zeros(length(signal),1);

[Nsamp,Nchannel] = size(signal);

t_add = zeros(Nsamp,1);

for n = 1:Nchannel % For each channel
    b = mean(signal(:,n)); 
    % Signal power for ECG is the poewr of the peak from the maximum to the
    % minimum point in the window around a peak for first 300 peaks and
    % taking average
    for i = 1:300 
        area = signal(ann(i+10)-window/2:ann(i+10)+window/2,n);
        pk_pk(i) = max(area)-min(area);
    end
    range = min(maxk(pk_pk,15));
    rang_e = max(mink(pk_pk,15));
    S = (pk_pk(pk_pk<range));
    S = mean(S>rang_e);

    % Take power of noise for first 300 samples
    for i = 1:300
        temp = noise(i*Fs+1:(i+1)*Fs+1,n);
        avg = mean(temp);
        dif = sqrt((avg-temp).^2);
    end
    range = min(maxk(dif,15));
    rang_e = max(mink(dif,15));
    RMS_noise = dif(dif<range);
    RMS_noise = mean(RMS_noise>rang_e);
    N = RMS_noise^2;

    % Calculate noise needed and add to signal
    SNR = SNR/10;
    reqd = 10^SNR;
    a = (S/N)/reqd;
    a = sqrt(a);
    new = signal(:,n) + a*noise(:,n);

    % Move the signal to shift back to same energy as before
    b = mean(new)-b;
    out(:,n) = new - b;

end

% Create file and Time Vector into wfdb format
tm = 0:1/Fs:length(signal)/Fs-1/Fs;
mat2wfdb(out,record_name,Fs);

end

