% Function:
%   Compute bandwidth with signal energy criteria
% Inputs:
%   Xn   -   input signal
%   fs   -   sample rate
%   spec -   signal energy specification
% Outputs:
%   BW   -   bandwidth
% Author: Jungwei Huang, NSYSU

function [BW] = bandwidth(Xn, fs, spec)

t = length(Xn)/fs; 
Xn = Xn - mean(Xn);
Xn_fft=fft(Xn, t*fs);
[theta, rho]=cart2pol(real(Xn_fft),imag(Xn_fft));

energy_ratio = floor(100* cumsum(rho(2:end/2+1).^2)/ sum(rho(2:end/2+1).^2));
BW = find(energy_ratio == spec, 1, "first")/t
