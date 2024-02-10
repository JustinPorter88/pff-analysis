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