clear all
close all
clc

[filename, pathname] = uigetfile('*.mat','Select LMS mat file(s)','MultiSelect','on');

if iscell(filename)
  fprintf(['' num2str(length(filename)) ' file(s) selected.\n'])
  for cntr = 1:length(filename)
    load([pathname,filename{cntr}]) ;
    [x,t,dt,fs,f]=ReadLMS(Signal_0,Signal_1) ;
    save([pathname 'Processed_' filename{cntr}],'x','t','dt','fs','f')
    fprintf('.')
    if cntr/20 == floor(cntr/20)
      fprintf('\n')
    end
  end
  fprintf('\nAll LMS Data were saved in Processed .m Files\n');
else
  load([pathname,filename]) ;
  [x,t,dt,fs,f]=ReadLMS(Signal_0,Signal_1) ;
  save([pathname 'Processed_' filename],'x','t','dt','fs','f')
end

