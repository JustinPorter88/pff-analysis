function [yt2, t2] = filterData(t,data,params)
%% Double Reverse Multimodal Decomposition (DRMD)
%    The code first expands the original data by system reconstruction and prediction.
%    Then Double Reverse Filter (DRF) is applied to the expanded signal.
%    Signal expansion function: 'ERA_SISO_decay.m'
%    DRF function: 'RevRevFilt.m'
%

%% Reverse, Filter (first time), and Reverse data
[data,~,t] = RevForFilt(data,t,params.freq,params.order,params.band);

%% Data Expansion (including identification and prediction)
PointIden_0=1; % default setting is 1
PointIden=1000; % number of data points that are used to identify
PointPredict=6000; % number of data points that will be predicted

data_0=data(PointIden_0:PointIden_0+PointIden); % data that are chosen to identify
t_0=t(PointIden_0:PointIden_0+PointIden); % time corresponding to the data chosen

nSRM=1; %% the order of the reconstructed system
ConstOrder=[1]; %#ok<NBRAK> % the order that needs to be constructed
dt = t(2)-t(1);
fs = 1/(t(2)-t(1));

[~,~,~,RespPred]=ERA_SISO_decay(data_0,nSRM,ConstOrder,fs,PointPredict);

RespPred=RespPred(end:-1:1);
TimePred=t_0(1)-dt:-dt:t_0(1)-dt-(length(RespPred)-1)*dt; % time corresponding to the predicted response
TimePred=TimePred(end:-1:1)';

data_ex=[RespPred;data(PointIden_0:end)]; % expanded data, including the original response and predicted one
t_ex=[TimePred;t(PointIden_0:end)]; % time corresponding to the expanded data

%% Forward Filter
[~,yt1,t_ex] = RevForFilt(data_ex,t_ex,params.freq,params.order,params.band);
temp1=find(t_ex>=0,1);
yt2=yt1(temp1:end);
t2=t_ex(temp1:end);
end