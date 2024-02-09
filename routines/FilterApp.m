function FilterApp(t,x) 
idx = find(t>0) ;
%% Calculate data to display
[ ~,X_abs,Fz1 ] = FFT_simp(x(idx(1):end),1/(t(2)-t(1)));

%% Initialize uifigure
uifig = uifigure('Name','Select Frequency Values');
            uifig.Color = [0.94 0.94 0.94];
            uifig.Position = [300 400 800 500];
uiax = uiaxes(uifig) ;
            xlabel(uiax, 'Frequency, Hz')
            ylabel(uiax, 'Response Amplitude')
            uiax.FontSize = 12;
            uiax.Box = 'on';
            uiax.XGrid = 'on';
            uiax.YGrid = 'on';
            uiax.Position = [25 200 750 275];
semilogy(uiax,Fz1,X_abs)
hold(uiax,'on')
uiax.XLim = [Fz1(1) ceil(Fz1(ceil(end/2)))];
            
%% Frequency and Bandwidth
            % Create Peak Frequency Slider Label
            FreqL = uilabel(uifig);
            FreqL.HorizontalAlignment = 'right';
            FreqL.Position = [10 170 100 15]; %[46 123 94 15];
            FreqL.Text = 'Peak Frequency';
            % Create Peak Frequency Slider
            Freq = uislider(uifig);
            Freq.Position = [120 175 515 3];
            Freq.Limits = [Fz1(1) Fz1(floor(end/2))];
            Freq.Value = floor(Fz1(floor(end/4))) ;
            % Create Peak Frequency Numeric Text Box
            FreqTB = uieditfield(uifig, 'numeric');
            FreqTB.ValueDisplayFormat = '%9.2f';
            FreqTB.Position = [675 160 100 22];
            FreqTB.Limits = Freq.Limits;
            FreqTB.Value = Freq.Value;
            
            % Create Bandwidth Slider Label
            BandL = uilabel(uifig);
            BandL.HorizontalAlignment = 'right';
            BandL.Position = [10 120 100 15];
            BandL.Text = 'Half Bandwidth';
            % Create Bandwidth Slider
            Band = uislider(uifig);
            Band.Position = [120 125 515 3];
            Band.Limits = [0 round(Fz1(end)/8)];
            Band.Value = 10 ;
            % Create Peak Frequency Numeric Text Box
            BandTB = uieditfield(uifig, 'numeric');
            BandTB.ValueDisplayFormat = '%9.2f';
            BandTB.Position = [675 110 100 22];
            BandTB.Limits = Band.Limits;
            BandTB.Value = Band.Value;

            % Create event handles
            FreqL = plot(uiax,Freq.Value*ones(100,1),linspace(uiax.YLim(1),uiax.YLim(2),100),'k') ;
            BandLU = plot(uiax,(Freq.Value+Band.Value)*ones(100,1),linspace(uiax.YLim(1),uiax.YLim(2),100),'--k') ;
            BandLL = plot(uiax,(Freq.Value-Band.Value)*ones(100,1),linspace(uiax.YLim(1),uiax.YLim(2),100),'--k') ;
            Freq.ValueChangedFcn = @(Freq,event) ChangeFreq(Freq,FreqL,Band,BandLU,BandLL,FreqTB) ;
            Freq.ValueChangingFcn = @(Freq,event) ChangeFreq(event,FreqL,Band,BandLU,BandLL,FreqTB) ;
            FreqTB.ValueChangedFcn = @(FreqTB,event) ChangeFreq(event,FreqL,Band,BandLU,BandLL,Freq) ;
            Band.ValueChangedFcn = @(Band,event) ChangeBand(Band,Freq,BandLU,BandLL,BandTB) ;
            Band.ValueChangingFcn = @(Band,event) ChangeBand(event,Freq,BandLU,BandLL,BandTB) ;
            BandTB.ValueChangedFcn = @(BandTB,event) ChangeBand(event,Freq,BandLU,BandLL,Band) ;

            
%% Order selection
            % Create Order Drop Down Label
            OrderDropDownLabel = uilabel(uifig);
            OrderDropDownLabel.HorizontalAlignment = 'right';
            OrderDropDownLabel.Position = [10 60 100 15];
            OrderDropDownLabel.Text = 'Filter Order';

            % Create Order Drop Down Selection
            OrderDD = uidropdown(uifig);
            OrderDD.Items = {'2', '3', '4', '5', '6', '7'};
            OrderDD.Position = [115 55 75 22];
            OrderDD.Value = '3';
            

%% Time cutoff selection
            % Create Cutoff Time Label
            CutoffLabel = uilabel(uifig);
            CutoffLabel.HorizontalAlignment = 'right';
            CutoffLabel.Position = [190 60 100 15];
            CutoffLabel.Text = 'Cutoff Time';
            % Create Cutoff Time Numeric Text Box
            CutoffTB = uieditfield(uifig, 'numeric');
            CutoffTB.ValueDisplayFormat = '%9.2f';
            CutoffTB.Position = [300 55 100 22];
            CutoffTB.Limits = [0 max(t)];
            CutoffTB.Value = floor(max(t)*10)/10 ;
            % Create Time Series Button
            CutoffButton = uibutton(uifig) ;
            CutoffButton.BackgroundColor = [0 0.451 0.7412];
            CutoffButton.FontWeight = 'bold';
            CutoffButton.FontColor = [1 1 1];
            CutoffButton.Position = [260 25 125 25];
            CutoffButton.Text = 'View Time Series';

            % Create exit handles
            CutoffButton.ButtonPushedFcn = @(CutoffButton,event) plotTS(t,x) ;            

            
%% Check Filter Results
            % Create Filter Check Button
            CheckButton = uibutton(uifig) ;
            CheckButton.BackgroundColor = [0 0.451 0.7412];
            CheckButton.FontWeight = 'bold';
            CheckButton.FontColor = [1 1 1];
            CheckButton.Position = [440 25 150 50];
            CheckButton.Text = 'Check Filter Results';

            % Create exit handles
            CheckButton.ButtonPushedFcn = @(CheckButton,event) FiltCheck(Freq,Band,OrderDD,CutoffTB,t,x) ; 


            
%% Exit button            
            ConfirmFrequenciesButton = uibutton(uifig) ;
            ConfirmFrequenciesButton.BackgroundColor = [0 0.451 0.7412];
            ConfirmFrequenciesButton.FontWeight = 'bold';
            ConfirmFrequenciesButton.FontColor = [1 1 1];
            ConfirmFrequenciesButton.Position = [635 25 140 50];
            ConfirmFrequenciesButton.Text = 'Confirm Parameters';

            % Create exit handles
            ConfirmFrequenciesButton.ButtonPushedFcn = @(ConfirmFrequenciesButton,event) FiltOut(Freq,Band,OrderDD,CutoffTB,uifig) ;
            uifig.CloseRequestFcn = @(uifig, event) FiltOut(Freq,Band,OrderDD,CutoffTB,uifig) ;
            

%% Command window and pause
fprintf('\nWaiting for user input through the graphical interface...\n')
waitfor(uifig) ;
end

%% Functions
function plotTS(t,x)
% Plot the time series so that the user can decide on a cutoff time
  figure
  plot(t,x)
  xlabel('Time')
  ylabel('Response')
  box on
end

%%
function ChangeBand(thing,Freq,line1,line2,other)
% Change the bandwidth display
  line1.XData = (thing.Value+Freq.Value)*ones(100,1);
  line2.XData = (-thing.Value+Freq.Value)*ones(100,1);
  other.Value = thing.Value ;
end

%%
function ChangeFreq(thing,line,Band,BandLU,BandLL,other)
% Change the frequency display
  line.XData = thing.Value*ones(100,1) ;
  other.Value = thing.Value ;
  ChangeBand(Band,thing,BandLU,BandLL) ;
end

%%
function FiltCheck(Freq,Band,OrderDD,CutoffTB,t,x)
% Test the filter parameters on the provided data
% Copy over parameters
params.freq = Freq.Value;
params.band = 2*ceil(Band.Value);
params.order = str2double(OrderDD.Value) ;
params.tlast = (CutoffTB.Value) ;

% Prep the data
idx = find(abs(diff(x)) == max(abs(diff(x))));
if idx > 200
  x = x - mean(x(1:round(idx/2))) ;
else
  x = x - mean(x(end-100:end)) ;
end
idx = find(abs(x)>0.5*max(abs(x))) ;
t = t - t(idx(1)-1) ;
idx = find(t>0,1,'first') ;
t = t(idx:end);
x = x(idx:end);

% %% Reverse and Filter data
% [x,~,t] = RevForFilt(x,t,params.freq,params.order,params.band);
%% Expand and Filter data
[x,t] = filterData(t,x,params) ;

% Cut data length
idx = find(t>params.tlast,1,'first') ;
t = t(1:idx) ;
x = x(1:idx) ;

% Extract the data
[Freq,Damp,Amp,Time]=FreqDampAmpExtraction(x,t);

% Plot the data
figure
subplot(131)
semilogx(Amp,Freq,'.');
xlabel('Amplitude');
ylabel('Frequency');
subplot(132)
semilogx(Amp,Damp*100,'.');
xlabel('Amplitude');
ylabel('Damping Ratio, %');
subplot(133)
plot(Time,Amp,'.');
xlabel('Time');
ylabel('Amplitude');
set(gcf,'Position',[325 450 1200 420])
end

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
ConstOrder=[1]; % the order that needs to be constructed
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


%%
function FiltOut(Freq,Band,OrderDD,CutoffTB,uifig)
% Export everything
global PFFparams

PFFparams.freq = Freq.Value;
PFFparams.band = 2*ceil(Band.Value);
PFFparams.order = str2double(OrderDD.Value) ;
PFFparams.tlast = (CutoffTB.Value) ;
selection = questdlg( ...
  'Exit filter parameter app and confirm chosen parameters?',...
    'Confirmation','Yes','No','Yes');
switch selection
  case 'Yes'
    delete(uifig)
  case 'No'
    return
end
end