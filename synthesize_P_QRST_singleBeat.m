% -------------------------------------------------------
%
%    synthesize_P_QRST_singleBeat.m  - This script synthesizes a simulated
%    P wave and a QRST segment to a full hearbeat
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


%% plotting
close all; 
leadID = 2; 

figure; 
hold on; 
plot(1:1:pwd_query+1, ecg_beat(leadID, 1:1:pwd_query+1), 'Linewidth', 2);
hold on; plot(pwd_query+1:1:pwd_query+poff_qon_int+1, ecg_beat(leadID,pwd_query+1:pwd_query+poff_qon_int+1), 'Linewidth', 2);
hold on; plot(pwd_query+poff_qon_int+1:1:pwd_query+poff_qon_int+qt_orig+1, ecg_beat(leadID,pwd_query+poff_qon_int+1:pwd_query+poff_qon_int+qt_orig+1), 'Linewidth', 2);
hold on; plot(pwd_query+poff_qon_int+1+qt_orig:1:rr_int, ecg_beat(leadID,pwd_query+poff_qon_int+1+qt_orig:rr_int), 'Linewidth', 2)


legend('p wave', 'pq segment', 'qrst segment', 'tp interval');
xlabel('time in ms');
ylabel('voltage in mV');
title('synthesized ECG from simualted P wave and QRST complex');
set(gca, 'FontSize', 20, 'fontname', 'times')
set(gcf, 'Position', [744   533   937   516]);