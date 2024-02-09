%% Use the PFF algorithm to analyze user selected data sets.
%  This code is based on the work by Mengshi Jin to develop the PFF
%  algorithm, and the work by Wei Chen and Mengshi Jin to develop the DRMD
%  filter algorithm.
%    Mengshi Jin               Wei Chen 
%    msjin@tongji.edu.cn       meshiawei@tongji.edu.cn 
%
%  GUI and user interaction coded by M. Brake (brake@rice.edu)
%  Revision date: 2020/12/07 - Minor changes, including handling transposed
%                              data sets
%  Revision date: 2020/05/27 - Incorporated the DRMD filter and 95%  
%                              confidence intervals
%  Revision date: 2020/05/14 - First interactive version of the code
%                              published

%#ok<*SAGROW>
%#ok<*CLALL>

clear all 
% close all
clc

addpath('./routines/')

%% Have the user specify the data sets of interest
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
waitfor(msgbox('\fontsize{11}Please select the data files to process. These should be *.mat files that contain Nx1 vectors with both a time record and data series.',...
    'Input Data Instructions','help',CreateStruct));
  
[filename, pathname] = uigetfile('*.mat','Select Data Files','MultiSelect','on') ;

if iscell(filename)
  load([pathname,filename{1}]) ;
else
  load([pathname,filename]) ;
end

%% Have the user specify the variables to be analyzed
[varout,varoutnames] = uigetvar({'Please select the time record.', ...
  'Please select the corresponding data signal'},'Introduction', ...
  'This UI requires that all data be stored in vectors. If you need to select a nested variable within a structure, or a particular column of a vector, please modify the program or data file.', ...
  'InputDimensions',[1 1]);

%% Check to see if a configuration file exists, otherwise use user input
waitfor(msgbox('\fontsize{11}Please select a .mat file with the filter parameters. This data file should contain a structure named params with children freq, band, order, and tlast that specify the center frequency, bandwidth, filter order, and length of time signal to analyze. See the source code for more details. If you do not already have a data file, select cancel on the next window, and a GUI will be opened to assist you in creating a set of filter specifications. After creating a filter, see the command window for information about where the filter has been saved for future use.',...
    'Filter Specification Instructions','help',CreateStruct));

[cfilename, cpathname] = uigetfile('*.mat','Select filter configuration file or cancel to create one') ;

if isequal(cfilename,0)

  fprintf('\nGenerating a new filter configuration file.')
%% Have the user specify the filter parameters
  global PFFparams %#ok<TLEV>

  [PFFtime,PFFdata] = prepData(varout{1},varout{2});
  FilterApp(PFFtime,PFFdata) ;
  params = PFFparams ;
  clear PFFtraw PFFxraw PFFparams
  fprintf('\nThe filter parameters are:\n')
  fprintf([' Center frequency: ' num2str(params.freq) '\n'])
  fprintf([' Bandwidth:        ' num2str(params.band) '\n'])
  fprintf([' Filter order:     ' num2str(params.order) '\n'])
  fprintf([' Cutoff time:      ' num2str(params.tlast) '\n\n'])

  configFile = ['CustomPFFConfiguration_' datestr(now,'yyyy_mmmm_dd')] ;
  fprintf(['\nFilter configuration file saved in the data directory as:\n' configFile '.m\n\n'])
  save([pathname configFile],'params')
else
%% Load existing configuration file
  load([cpathname,cfilename],'params') ;
  fprintf('\nParamters successfully loaded:\n')
  fprintf([' Center frequency: ' num2str(params.freq) '\n'])
  fprintf([' Bandwidth:        ' num2str(params.band) '\n'])
  fprintf([' Filter order:     ' num2str(params.order) '\n'])
  fprintf([' Cutoff time:      ' num2str(params.tlast) '\n\n'])
end



%% Analyze the data using the PFF algorithm
%% For multiple files, perform the uncertainty analysis
if iscell(filename)
  % Analyze and save in big matrices
  for cntr = 1:length(filename)
    load([pathname,filename{cntr}]) ;
    eval(['PFFtimevec = ' varoutnames{1} ';']);
    eval(['PFFxvec = ' varoutnames{2} ';']);
    [Freq,Damp,Amp,Time] = PFF_Analyze(PFFtimevec,PFFxvec,params) ;
    Freqs{cntr} = Freq ;
    Amps{cntr} = Amp ;
    Damps{cntr} = Damp ;
    Times{cntr} = Time ;
  end

  % Calculate the average and 95% confidence interval
  [AMP_avg,FRE_avg,FRE_upper,FRE_lower,DAM_avg,DAM_upper,DAM_lower]=AvgCiC(Freqs,Amps,Damps);

  % Plot results
  color=[{'r'},{'k'},{'b'}];
  linestyle=[{'-'},{'--'},{':'}];
  linewidth=[{0.8},{2},{2}];

  figure
  hold on;grid on;box on;
  for j=1:length(Amps)
      str=strcat('plot(Amps{j},Freqs{j},',mat2str(color{1}),')');
      eval(str);
  end
  h1=area(AMP_avg,FRE_upper,'FaceColor','k','EdgeColor','w');
  h2=area(AMP_avg,FRE_lower,'FaceColor','w','EdgeColor','w');
  set(h1,'FaceAlpha',0.2);
  set(h2,'FaceAlpha',1);
  str=strcat('plot(AMP_avg,FRE_avg,''linestyle'',',mat2str(linestyle{2}),',''color'',',mat2str(color{2}),...
      ',''linewidth'',',mat2str(linewidth{2}),')');
  eval(str);
  str=strcat('plot(AMP_avg,FRE_upper,''linestyle'',',mat2str(linestyle{3}),',''color'',',mat2str(color{3}),...
      ',''linewidth'',',mat2str(linewidth{3}),')');
  eval(str);
  str=strcat('plot(AMP_avg,FRE_lower,''linestyle'',',mat2str(linestyle{3}),',''color'',',mat2str(color{3}),...
      ',''linewidth'',',mat2str(linewidth{3}),')');
  eval(str);
  set(gca,'xscale','log')
  set(gca,'layer','top');
  xlabel('Amplitude')
  ylabel('Frequency')

    figure
  hold on;grid on;box on;

  for j=1:length(Amps)
      str=strcat('plot(Amps{j},Damps{j}*100,',mat2str(color{1}),')');
      eval(str);
  end
  h1=area(AMP_avg,DAM_upper*100,'FaceColor','k','EdgeColor','w');
  h2=area(AMP_avg,DAM_lower*100,'FaceColor','w','EdgeColor','w');
  set(h1,'FaceAlpha',0.2);
  set(h2,'FaceAlpha',1);
  str=strcat('plot(AMP_avg,DAM_avg*100,''linestyle'',',mat2str(linestyle{2}),',''color'',',mat2str(color{2}),...
      ',''linewidth'',',mat2str(linewidth{2}),')');
  eval(str);
  str=strcat('plot(AMP_avg,DAM_upper*100,''linestyle'',',mat2str(linestyle{3}),',''color'',',mat2str(color{3}),...
      ',''linewidth'',',mat2str(linewidth{3}),')');
  eval(str);
  str=strcat('plot(AMP_avg,DAM_lower*100,''linestyle'',',mat2str(linestyle{3}),',''color'',',mat2str(color{3}),...
      ',''linewidth'',',mat2str(linewidth{3}),')');
  eval(str);
  set(gca,'xscale','log')
  set(gca,'layer','top');
  xlabel('Amplitude')
  ylabel('Damping Ratio, %')
else
  
%% For a single file, extract the data and plot  
  PFFtimevec = varout{1};
  PFFxvec = varout{2} ;
  [Freq,Damp,Amp,Time] = PFF_Analyze(PFFtimevec,PFFxvec,params) ;

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





function [Freq,Damp,Amp,Time] = PFF_Analyze(t_raw,x_raw,params) 
%% Prep data
[t,x] = prepData(t_raw,x_raw);
idx = find(t>0,1,'first') ;
if idx > 1
  t = t(idx-1:end);
  x = x(idx-1:end);
end

% %% Reverse and Filter data
% [x,~,t] = RevForFilt(x,t,params.freq,params.order,params.band);
%% Expand and Filter data using the DRMD method
[x,t] = filterData(t,x,params) ;

%% Cut data length
if ~isnan(params.tlast)
  idx = find(t>params.tlast,1,'first') ;
else
  if floor(t(end)*10) >= 15
    idx = find(t>floor(t(end)*10)/10,1,'first') ;
  else
    idx = find(t>floor(t(end)*100)/100,1,'first') ;    
  end  
end
t = t(1:idx) ;
x = x(1:idx) ;

%% Extract properties
[Freq,Damp,Amp,Time]=FreqDampAmpExtraction(x,t);

clear x t idx
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


function [PFFtime, PFFdata] = prepData(t,x)
sz = size(t);
if sz(1) < sz(2)
  t = t';
  x = x';
end

%% Orient the data to start at t = 0 for the start of the impact response.
idx = find(abs(diff(x)) == max(abs(diff(x))));
if idx > 200
  PFFdata = x - mean(x(1:round(idx/2))) ;
else
  PFFdata = x - mean(x(end-100:end)) ;
end
idx = find(abs(PFFdata)>0.5*max(abs(PFFdata))) ;
PFFtime = t - t(idx(1)-1) ;
end