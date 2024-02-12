function [Amp, Freq, Damp, Time] = pff_function(t, x, params)
%
% Function that can be called to apply Peak Finding and Fitting analysis
% without using the GUI. 
% 
% Inputs:
%   t - time vector to analyze
%   x - response vector to analyze 
%   params - structure with the parameters
%       .freq - center frequency for the filter (Hz)
%       .band - (full) bandwidth of filter in units of frequency (Hz)
%       .order - order of butterworth filter to use
%       .tlast - final time to cut off in analysis (only signal before this
%                time is used). This must be less than the total time or
%                the function will hit an error.
%
% Outputs:
%   Amp - amplitude of response vector (after filtering) at discrete points
%   Freq - calculated frequency at same discrete points
%   Damp - calculated damping at same discrete points
%   Time - times of the calculated points. 
%
% Notes:
%   1. This script can be verified with the script
%   tests/test_pff_function.m
%   2. This script (should) automatically add the path of the routines
%   folder to your MATLAB path so that it can execute PFF.

    
    % Add the path to the dependencies
    [function_path, ~, ~] = fileparts(which('pff_function'));
    dependencies_path = fullfile(function_path, 'routines');
    addpath(dependencies_path)
    addpath(fullfile(dependencies_path, 'copied_from_gui'))

    % Do PFF Calculation
    [Freq,Damp,Amp,Time] = PFF_Analyze(t,x,params) ;

end