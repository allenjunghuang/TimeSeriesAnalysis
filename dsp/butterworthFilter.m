% Function: 
%   None phase shift Butterworth filtering
% Inputs:
%   Xn   -   input signal
%   fs   -   sample rate
% Outputs:
%   Yn   -   output signal
%   fc   -   cutoff frequency
% Reference: 
% [1] B.P. Lathi. Signal Processing and Linear System, Oxford
% Author: Jungwei Huang, NSYSU

function [Yn, fc, num, den] = buterworthIIR (Xn, pass_hz, stop_hz, fs)

if (size(Xn,1) > size(Xn,2)) %correct data dimension to row vector
    S1 = Xn';
else
    S1 = Xn;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all; clc; close all;
% T0=1; %the signal interval we interested
% fs=512; %sampling frequency
% T=1/fs; %sampling interval
% N0=T0/T; % number of samples in T0
% t=0:T:T*(N0-1); 
% S1=cos(10*2*pi*t)+cos(50*2*pi*t);
% pass_hz=25;%[25 30];
% stop_hz=30;%[22 33];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Butterworth Filter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_pass = (pass_hz)./(0.5*fs); 
w_stop = (stop_hz)./(0.5*fs); 
rpl_band = 2; %passband ripple
atn_stop = 5; %stopband attenuation
[n, cutoff_hz] = buttord(w_pass, w_stop, rpl_band, atn_stop); %for lowpass and bandpass filter. note the unit is pi.
% [z,p,k]=tf2zp(num,den);
if length(w_pass) ~= length(w_stop) %% check Butterworth filter settings
    fprintf('Butterworth filter setting error');
else
    fc = cutoff_hz*(0.5*fs); 
    if length(w_pass) > 1 && length(w_stop) > 1 
        if w_stop(2) > w_pass(2) > w_pass(1) > w_stop(1)
            [num, den] = butter(n, cutoff_hz);
            fprintf('Butterworth bandpass filter >> order: %d,cutoff frequency: %6.2f Hz %6.2f Hz \n',n,fc(1),fc(2));
        elseif w_pass(2) > w_stop(2) > w_stop(1) > w_pass(1)
            [num, den] = butter(n, cutoff_hz,'stop');
            fprintf('Butterworth bandstop filter >> order: %d,cutoff frequency: %6.2f Hz %6.2f Hz \n',n,fc(1),fc(2));
        end
    elseif length(w_pass) == 1 && length(w_stop) == 1        
        if w_stop(1) > w_pass(1)
            [num, den] = butter(n,cutoff_hz,'low');
            fprintf('Butterworth lowpass filter >> order: %d,cutoff frequency: %6.2f Hz \n',n,fc);
        elseif w_pass(1)>w_stop(1)
            [num, den] = butter(n,cutoff_hz,'high');
            fprintf('Butterworth highpass filter >> order: %d,cutoff frequency: %6.2f Hz \n',n,fc);
        end        
    else
       fprintf('Butterworth filter setting error');
    end
end

% [h, w]=freqz(num,den,fs); %buterworth filter frequency response
% figure(123) %buterworth filter spectrum
% subplot(2,1,1); plot(w/pi*fs/2, (abs(h)));
% xlabel('Hz'); ylabel('|D_n|');
% title('Butterworth filter spectrum'); grid on;
% subplot(2,1,2); plot(w/pi*fs/2, unwrap(angle(h))*180/pi); 
% grid on;
% xlabel('Hz'); ylabel('\angle D_n (deg)');

% S1_fft=fft(S1);
% pDn_S1=S1_fft./N0; 
% Dn_S1=pDn_S1(1:(N0/2));
% [Dnangle_S1, Dnmag_S1]=cart2pol(real(Dn_S1), imag(Dn_S1));
% w_seq=0:(N0/2-1);

S2=fliplr(S1); %flip matrix left to right
S3=filter(num,den, S2); %1st filter
S4=fliplr(S3);
Yn=filter(num,den, S4); %2ed filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier Spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% S4_fft=fft(S4); %1st filtered signal FFT
% pDn_temp=S4_fft./N0; 
% Dn_temp=pDn_temp(1:(N0/2));
% [Dnangle_temp,Dnmag_temp]=cart2pol(real(Dn_temp),imag(Dn_temp));
% 
% output_fft=fft(output_2); %2ed filtered signal FFT
% pDn_output=output_fft./N0; 
% Dn_output=pDn_output(1:(N0/2));
% [Dnangle_output,Dnmag_output]=cart2pol(real(Dn_output),imag(Dn_output));
% 
% figure(101)
% plot(t,S1,'r'); hold on;
% plot(t,S4,'g'); hold on;
% plot(t,Yn,'b');
% xlabel('t (sec)'); ylabel('f(t)');
% legend('Original: cos(20\pi t)+cos(100\pi t)','First Filtered','Second Filtered');
% 
% figure(102)
% subplot(2,1,1);
% plot(w_seq,Dnmag_S1);
% title('Original Spectrum');
% grid on; xlabel('Hz'); ylabel('|D_n|');
% subplot(2,1,2); plot(w_seq, unwrap(Dnangle_S1)*180/pi); grid on;
% xlabel('Hz'); ylabel('\angle D_n (deg)');
% 
% figure(103)
% subplot(2,1,1); 
% plot(w_seq,Dnmag_temp);
% title('First Filtered Specturm');
% grid on; xlabel('Hz'); ylabel('|D_n|');
% subplot(2,1,2); plot(w_seq, unwrap(Dnangle_temp)/pi*180); grid on;
% xlabel('Hz'); ylabel('\angle D_n (deg)');
% 
% figure(104)
% subplot(2,1,1);
% plot(w_seq,Dnmag_output);
% title('Second Filtered Specturm');
% grid on; xlabel('Hz'); ylabel('|D_n|');
% subplot(2,1,2); plot(w_seq, unwrap(Dnangle_output)/pi*180); grid on;
% xlabel('Hz'); ylabel('\angle D_n (deg)');

return