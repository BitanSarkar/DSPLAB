clc;
clear all;
close all;

data = load('Exp03_PPG_25hz_75samples.csv');
x = data; window = 8;N = size(x,2);
sm = zeros(size(data));
fs = 25;

% Moving Average Filter
for i = 1:N
    s = 0;
    if i <= window
        s = sum(x(1:i));
    else
        s = sum(x(i-window:i));
    end
    sm(i) = s / window;
end

% DFT and DFT Matrix
DFT = fft(sm); DFT_coef = dftmtx(N);
[~,loc] = max(abs(DFT(3:end/2)));
mag = abs(DFT);


figure;
plot(mag);
mag = abs(DFT);xlabel('Frequency(in Hz)'); ylabel('Magnitude');axis tight;

fu = 0.5;
si = round( (fu*N)/ fs);

figure;
for i= si:N-si;
    mag(i)=0;
end
plot(mag); xlabel('Frequency(in Hz)'); ylabel('Magnitude');axis tight;
title('Magnitude Response'); axis tight; grid on;

re = ifft(mag);

figure;
plot(abs(re));

% Pulse Rate Estimation using DFT
display('Pulse Rate using DFT :');
display(60*(loc+1)* fs/N);

% Autocorrelation
sm = sm - mean(sm);
sm = sm/abs(max(sm));
acf = xcorr(sm);
[maxa,acf_loc]=max(acf);
acf = acf(acf_loc+1:end);
acf=acf/maxa;

%First Zero Crossing 
fzcp = 0;
for i=1:size(acf, 2)
    if acf(i+1) * acf(i) < 0
        fzcp = i + 1;
        break
    end
end

% Pulse Rate using Autocorrelation
[maxi, maxi_loc] = max(acf(fzcp:end));
maxi_loc = fzcp + maxi_loc;
display('Pulse Rate using ACF sequence :');
display(60*fs/(maxi_loc-1));
