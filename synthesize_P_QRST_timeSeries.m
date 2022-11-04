% -------------------------------------------------------
%
%    synthesize_P_QRST_timeSeries.m  - This script synthesizes a simulated
%    P wave and a QRST segment to a full hearbeat and extends it to a 10s
%    time series
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
% Revision history:
%
%
%% stitch load relevant data and add paths
close all; clear all; 
% add path to ECGdeli (to be downloaded here: https://github.com/KIT-IBT/ECGdeli)
addpath(genpath('../ECGdeli'));

% load multivariate normal distributions for amplitudes and intervals
load('multivariate_distributions/MVD.mat')

% load simulated p waves and qrs complexes (unfiltered!)
load('example_signals/p1_sim_unfiltered.mat');
load('example_signals/qrst1_sim_unfiltered.mat');

% load noise vector
flag_noise = false;
if flag_noise == true
    load('noise_data/DATA_noises_real.mat');
end

%% scale QRST segments according to p wave amplitude in II (main heart axis!)
query_pamp = max(p(2,:)); 

% draw random samples from amplitude distributions
rnd_amps = mvnrnd(GMModel_amp.mu,GMModel_amp.Sigma,10000000);

% find all values in clinical amplitude distribution that have a similar p
% wave amplitude than the simulated signal
val = find(abs(rnd_amps(:,1)-query_pamp)<5e-5);
R_constraint_amp = rnd_amps(val,2);

% use the mean R peak value fitting to the P peak amplitude - alternative
% approach: draw random sample from 1 dimensional R_constraint_amp distribution
val_sample_a = mean(R_constraint_amp);
scaling_qrs = val_sample_a/max(qrst(2,:));

% apply scaling.
qrst_scaled = qrst.*scaling_qrs;

%% stitch scaled signals together

% apply ecg deli to find pwd and qt intervals

% two methods to stitch p and qrst together: sigmoid function or straight
% line - line is used in the example below
sigm = @(x) 1./(1+exp(-a*x));
line = @(x, m) x.*m;

% apply ECGdeli and find annotations
ecg_tmp = [];
for leadID = 1:12
    m_tmp1 = (qrst_scaled(leadID,1)-p(leadID,end))/(100);
    m_tmp2 = (p(leadID,1)-qrst_scaled(leadID,end))/(200);
    ecg_tmp(leadID,:) = [p(leadID,:), line(1:1:100,m_tmp1), qrst_scaled(leadID,:)-p(leadID,end), line(1:1:200,m_tmp2)];
end
ecg_tmp_rep = repmat(ecg_tmp, 1,20);

[FPT_MultiChannel,~]=Annotate_ECG_Multi(ecg_tmp_rep',1000);

% find p wave duration -> use second beat, since the first one is often
% prone to errors.
pwd_query = FPT_MultiChannel(2,3)-FPT_MultiChannel(2,1);

% cut out the p wave found by ECGdeli.
p_crop = p(:,1:pwd_query);

% find qt duration
qt_query = FPT_MultiChannel(2,12)-FPT_MultiChannel(2,4); 
qt_orig = size(qrst_scaled,2);

% generate random samples from the intervals' multivariate distribution

rnd_int = mvnrnd(GMModel_int.mu,GMModel_int.Sigma,10000000);

% find all random numbers close to the simulated query values
rnd_qt = find(abs(rnd_int(:,3)-qt_query)<3);
rnd_pwd = find(abs(rnd_int(:,2)-pwd_query)<3);
val=intersect(rnd_qt,rnd_pwd);

R_constraint = rnd_int(val,[1,4]); 

GMModel_constraint_int = fitgmdist(R_constraint,1);

poff_qon_int = round(GMModel_constraint_int.mu(1,1));
rr_int = round(GMModel_constraint_int.mu(1,2));


ecg_beat = qrst_scaled(:,end).*ones(12,rr_int);

for leadID = 1:12
    m = (qrst_scaled(leadID,1)-p_crop(leadID,end))/(poff_qon_int);
    pq_sig = p_crop(leadID,end) + line(1:1:poff_qon_int, m);
    ecg_beat(leadID,1:pwd_query) = p_crop(leadID,:);
    ecg_beat(leadID,pwd_query+1:pwd_query+poff_qon_int) = pq_sig;
    ecg_beat(leadID,pwd_query+poff_qon_int+1:pwd_query+poff_qon_int+qt_orig) = qrst_scaled(leadID,:);
    
    % prepare for the next beat -> ensure that qrst part ends at the same
    % signal value than p beginning
    tp_dur = rr_int - (pwd_query+poff_qon_int+qt_orig); 
    m = (p_crop(leadID,1)-qrst_scaled(leadID,end))/(tp_dur);
    tp_sig = qrst_scaled(leadID,end) + line(1:1:tp_dur, m);
    
    ecg_beat(leadID,pwd_query+poff_qon_int+1+qt_orig:end) =  tp_sig;

end

%% generate RR time series
[tau] = 1000.*f_getRRseries(rr_int*0.001, 15);
qt_mod = nan(size(tau,1),1);

% find QT intervals for RR time series. 
for t_ID = 1:size(tau,1)
    % draw random samples from amplitude distributions
    rnd_int = mvnrnd(GMModel_int.mu,GMModel_int.Sigma,10000000);

    % find all random numbers close to the simulated query values
    rnd_rr = find(abs(rnd_int(:,4)-tau(t_ID,1))<3);
    rnd_pqi = find(abs(rnd_int(:,1)-poff_qon_int)<3);
    rnd_pwd = find(abs(rnd_int(:,2)-pwd_query)<3);
    val=intersect(rnd_rr,rnd_pwd);

    R_constraint = rnd_int(val,[3]); 

    GMModel_constraint_int = fitgmdist(R_constraint,1);

    qt_mod(t_ID,1) = round(GMModel_constraint_int.mu(1,1));
    tp_mod(t_ID,1) = tau(t_ID,1)-poff_qon_int-pwd_query-qt_mod(t_ID,1);
    
end


%% put time series everything together

% take the single beat and change QT segment as well as TP segment.

for t_ID = 1:size(tau,1)
    for leadID = 1:12
        qt_new = resample(qrst_scaled(leadID,:),qt_mod(t_ID,1),qt_query);
        tp_new  = tau(t_ID,1) - poff_qon_int-pwd_query-size(qt_new,2);
        m = (p_crop(leadID,1)-qrst_scaled(leadID,end))/(tp_new);
        tp_sig = qrst_scaled(leadID,end) + line(1:1:tp_new, m);
    
        ecg_beats{t_ID,leadID} = [ecg_beat(leadID,1:pwd_query+poff_qon_int), qt_new, tp_sig];
    end
    
end
ecg_series = [];
for leadID = 1:12
    ecg_series(leadID,:) = cell2mat(horzcat(ecg_beats(:,leadID)'));
end

%% add noise and filter signals

if flag_noise == true
    start_noise = randi(size(mixture_of_noises,2)-size(ecg_series,2)-1);
    noise_amp = 0.03;
    ecg_series_noise = ecg_series+noise_amp.*mixture_of_noises(1:12,start_noise+1:start_noise+size(ecg_series,2)); 
else
    ecg_series_noise = ecg_series;
end
[ECGsignals_filtered] = f_filterECG(ecg_series_noise', 1000, 0.5, 150);

%% plotting

leadIDs = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
figure; hold all; 
dist = 2; 
for leadID = 1:12
    % synethsized
    sig = (ECGsignals_filtered(leadID,1:10000)-ECGsignals_filtered(leadID,1))-leadID*dist;
    
    sig = sig(1,1:10000); % only select the first 10 seconds
    plot(sig, 'k', 'LineWidth', 1);    
    txt = leadIDs{1,leadID};
    text(10,sig(1,10)+0.15,txt)
end
title('Synthesized 10s ECG')
xlabel('time in ms')
ylabel('voltage in mV')
ylim([-25 1]);
yticks();
grid minor
set(gcf, 'Position',  [331         264        2173        1054])

%% save ECG in json file format
f_writeJson(ECGsignals_filtered, 1000, 'example_signals/p1_qrst1_synthesized.json');
