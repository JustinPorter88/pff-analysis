function [freq,freq_t]=FreqCal(T_max,T_min)
%% Frequency calculation in the PFF algorithm using backward difference
%
% input:
% T_max: time instants of fitted local peaks
% T_min: time instants of fitted local valleys
%
% output:
% freq: instantaneous frequency
% freq_t: corresponding time of the instantaneous frequency

%%
T=[T_max;T_min];
T=sort(T);
kc=1;
for loop1=1:length(T)-1
    freq(kc,1)=1/(T(loop1+1)-T(loop1))/2;
    kc=kc+1;
end
freq_t=T(1:end-1);

end