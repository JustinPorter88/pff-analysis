function [Freq,Damp,Amp,Time]=FreqDampAmpExtraction(data,t)
%% Identification of instantaneous frequency, damping ratio, and amplitude using PFF algorithm

% note: the calculated frequency and damping are based on both the upper and lower amplitude
%
% input: (with the same length)
% data: data
% t: time
%
% output: (with the same length)
% Freq: raw instantaneous frequency (Hz)
% Damp: raw instantaneous damping ratio
% Amp: raw instantaneous amplitude
% Time: corresponding time

%%
[X_max,T_max,X_min,T_min,~,~]=PeakFinding(data,t);
[freq_temp,~]=FreqCal(T_max,T_min);
[A_temp,T_temp]=com_max_min(X_max,T_max,X_min,T_min);
A_temp=A_temp(1:end-1);
T_temp=T_temp(1:end-1);
[Damp,Freq,Amp,Time]=DampCal_c(freq_temp,A_temp,T_temp);

end