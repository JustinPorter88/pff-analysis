% This script is to be used to run a test to verify that the pff_function
% correctly matches the previous GUI results. 
%
% Data to be used:
%   ../data/Sample_LMS_Data.mat - data that is processed (x versus t)
%   ../data/Sample_LMS_Params.mat - parameters for PFF to be used
%   ../data/Sample_LMS_Results.mat - reference results produced using
%   PFF_Analysis GUI.

addpath('..')

%% Load Data

data = load('../data/Sample_LMS_Data.mat', 't', 'x');
params = load('../data/Sample_LMS_Params.mat', 'params').params;
reference = load('../data/Sample_LMS_Results.mat');

%% Apply PFF

[Amp, Freq, Damp, Time] = pff_function(data.t, data.x, params);

%% Calculate Errors

Amp_err = max(abs(reference.Amp - Amp));
Freq_err = max(abs(reference.Freq - Freq));
Damp_err = max(abs(reference.Damp - Damp));
Time_err = max(abs(reference.Time - Time));

Tot_err = Amp_err + Freq_err + Damp_err + Time_err;

if Tot_err == 0
    disp('Test passed.');
else
    assert(false, 'Test failed! - Investigate the sources of error.')
end
