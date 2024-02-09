function [yt]=RevRevFilt(data,t,fn,n_bp,bw_bp,ftype)
%% Double Reverse Filter
%
% input:
% data: orignal data to be filtered
% t: time
% fn: central frequency of the filter (Hz)
% n_bp: order of the filter
% bw_bp: bandwidth of the filter
% ftype: type of the filter, 'bandpass' or 'lowpass' (for low frequency)
%
% output:
% yt: double reverse filtered signal (may have data before t=0)

%%
fs = 1/(t(2)-t(1));
if strcmp(ftype,'bandpass')
    fc_bp = [fn-bw_bp/2,fn+bw_bp/2 ]/fs*2;
else if strcmp(ftype,'low')
        fc_bp = (fn+bw_bp)/fs*2;
    end
end

[b,a]=butter(n_bp,fc_bp,ftype);
ut=data(end:-1:1);
vt = filter(b,a,ut);
wt=vt(end:-1:1);
yt = filter(b,a,wt);

end

