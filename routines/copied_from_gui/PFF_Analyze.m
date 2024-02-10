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