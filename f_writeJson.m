% -------------------------------------------------------
%
%    f_writeJson - Write ECG data into json files
%
%    Ver. 1.0.0
%
%    Created:           Claudia Nagel (13.10.2022)
%    Last modified:     Claudia Nagel (13.10.2022)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
%
% function [] = f_writeJson(signal, Fs, filename)
%
% Inputs:
%       signal: ECG signal (LxT vector, L: leads, T: timesteps)
%       Fs: sampling rate
%       filename: filename the data should be stored in
%
% Outputs:
%
%
% Example Usage:
%       f_writeJson(ecg, 500, 'ecg_synthesized.json')
%
% Revision history:
%
%

function [] = f_writeJson(signal, Fs, filename)

% specify all fields in json struct
val.t = 0:1:size(signal,2)-1;
val.ecg.I = signal(1,:);
val.ecg.II = signal(2,:);
val.ecg.III = signal(3,:);
val.ecg.aVR = signal(4,:);
val.ecg.aVL = signal(5,:);
val.ecg.aVF = signal(6,:);
val.ecg.V1 = signal(7,:);
val.ecg.V2 = signal(8,:);
val.ecg.V3 = signal(9,:);
val.ecg.V4 = signal(10,:);
val.ecg.V5 = signal(11,:);
val.ecg.V6 = signal(12,:);
val.Fs = Fs; 
json_struct = jsonencode(val);

fid = fopen(filename,'wt');
fprintf(fid, json_struct);
fclose(fid);
end

