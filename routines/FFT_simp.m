function [ X,X_abs,Fz ] = FFT_simp( x,fs )
%% Inputs:
% Response vector x
% Sampling frequency fs
%
%% Outputs:
% Response and response magnitude: X and X_abs
% Frequencies: Fz

L=length(x);
N = 2^nextpow2(L); % FFTs require 2^n samples, where n is an integer
Y = fft(x,N)/L*2;  % Built in Matlab FFT command
Y(1) = Y(1)/2;     % From Y, calculate the single sided response spectrum X 
X=Y(1:N/2+1);      
X_abs=abs(X);
Fz = fs/2*linspace(0,1,N/2+1);
end

