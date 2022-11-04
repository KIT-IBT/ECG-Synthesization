% -------------------------------------------------------
%
%    f_filterECG.m  - This function applies high- and lowpass filters and
%    removes isoline offset from ECGs
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
% function [ECGsignals_filtered] = f_filterECG(ECGsignal, samplerate, hp_frq, lp_frq)
%
% Inputs:
%       signal: ecg signal to be filtered
%       samplerate: sample frequency in Hz of the signal 
%       lp_frq: cutoff frequency of the lowpass filter 
%       hp_frq: cutoff frequency of the highpass filter 

%
% Outputs:
%       ECGsignals_filtered: filtered ecg signal.
%
%
% Example Usage:
%       ecg_filtered = f_filterECG(ecg, 500, 0.5, 150)
%
% Revision history:
%
%

function [ECGsignals_filtered] = f_filterECG(ECGsignal, samplerate, hp_frq, lp_frq)

[ecg_filtered_frq]=ECG_High_Low_Filter(ECGsignal,samplerate,hp_frq,lp_frq);
[ECGsignals_filtered,~,~,~]=Isoline_Correction(ecg_filtered_frq);
%ECGsignals_filtered=Notch_Filter(ecg_filtered_frq,samplerate,50,1);
ECGsignals_filtered = ECGsignals_filtered';
end

