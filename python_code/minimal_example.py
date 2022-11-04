
'''
--------------------------------------------------------------


Created:                          Nadine Eichhorn (28.08.2022)
 
 
 
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
 
This is a minimal exapmle for the code generateECG  


---------------------------------------------------------------      
 
Input Parser:
       
    1. P Wave: give the path to the file for the P wave
    2. QRST complex: give the path to the QRST complex file
    3. Folder to save the ECGs: give the path to the ordner where the ECGS should be saved
 
Example:
    
    1. example_signals/p2_sim_unfiltered.mat
    2. example_signals/qrst1_sim_unfiltered.mat
    3. Volumes/.../Resuls
 
 
Output:
   
    - 3 ECGs with 15 heartbeats:
        1. unfilterd ECG
        2. ECG with noise
        3. filterd ECG
    - Plot of the P wave
    - Plot of the QRST complex
    - Plot of one ECG

 The present working directory must be cn617_ecgsnthesization

 Example: python3 python_code/minimal_bsp.py  example_signals/p2_sim_unfiltered.mat example_signals/qrst1_sim_unfiltered.mat /Volumes/koala/Users/ne725/Documents/Hek/Ergebnis
 
---------------------------------------------------------------------                                            )
'''

# import used modules
 
from scipy.io import loadmat
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy import signal
import json
import csv
import argparse
import random
 
#Initialisiere Parser
 
parser=argparse.ArgumentParser(description='Load P wave and QRST complex. Enter a save location')
parser.add_argument('P_Welle',type=str, help='Give the path to the folder for the p waves')
parser.add_argument('QRST_Welle', type=str, help='Give the path to the folder for the qrst complexes')
parser.add_argument('Speicherort',type=str, help='Give the path to the folder in which the ECGs should be saved')
 
args=parser.parse_args()
 
 
'''
1. Dataloader
 
'''
 
# load p wave
 
p_wave=loadmat(args.P_Welle)
ecg_temp_rep_filt=[[element for element in upper]for upper in p_wave['p']]
p=np.array(ecg_temp_rep_filt)*1.5
 
# plot p wave
 
plt.plot(p[0])
plt.title('P Welle Ableitung II')
plt.show()
p[6,:]=p[6,:]*0.6
p[7,:]=p[7,:]*0.75
 
# plot qrst complex
 
qrst_wave=loadmat(args.QRST_Welle)
qrst_list=[[element for element in upper]for upper in qrst_wave['qrst']]
qrst=np.array(qrst_list)
 
#plot qrst complex
 
plt.plot(qrst[1,:])
plt.title('QRST 2.Ableitung')
plt.show()
 
 
'''
 
2. Single Beat Repition:
   
    -> scale QRST segments according to p wave amplitude in II (main heart axis!)
    -> stich the qrst complex and the p wave togehter
    -> use the single beat ECG to find the p and qt wave duration
   
 
'''
query_pamp= np.max(p[1,:])
 
# create an Multivariate Gaussian Distribution for the amplitudes
mean_amp = np.array([0.0894,0.8264])
sigma_amp = np.array([[0.0022,-5.8595e-04],[-5.8595e-04, 0.1101]])
 
rnd_amp = np.random.multivariate_normal(mean_amp, sigma_amp, 10000000)  # draw random samples from amplitude distributions
 
#find all values in clinical amplitude distribution that have a similar  pwave amplitude than the simulated signal
 
 
erfolg=False
wert=5e-5
 
while erfolg==False:
   
    val= np.nonzero(abs(rnd_amp[:,0]-query_pamp)<wert)
 
    wert=wert*10
 
    if len(val[0])==0:
      
        erfolg=False
        
    else:
 
        erfolg=True
       
 
R_constraint_amp= rnd_amp[val,1]
# use the mean R peak value fitting to the P peak amplitude - alternative approach: draw random sample from 1 dimensional R_constraint_amp distribution
 
val_sample_a=np.mean(R_constraint_amp)
scaling_qrs=val_sample_a/np.max(qrst[1,:])
 
# apply scaling
qrst_scaled= qrst*scaling_qrs
 
# stitch scaled signals through two straight lines
ecg_temp=[0]*12
for i in range(0,12):
   
    m_temp1=(qrst_scaled[i,0]-p[i,-1])/100                             # slope for the straight line from the end of the p wave to the beginning of the q wave
    m_temp2=(p[i,0]-qrst_scaled[i,-1])/200                             # slope for the straight line from the end of the qrst complex to the beginning of the p wave
   
    gerade1=np.arange(1,101)*m_temp1                                   # straight line to connect the end of the p wave to the beginning of the q wave
    gerade2= np.arange(1,201)*m_temp2                                  # straight line to connect the end of the ecg to the beginning of the new ecg
   
    qrst_lst=np.subtract(qrst_scaled[i,:],p[i,-1]).tolist()
    p_lst=p[i,:].tolist()
    p_lst.extend(gerade1)
    p_lst.extend(qrst_lst)
    p_lst.extend(gerade2)
   
    ecg_temp[i]=p_lst
      
 
# find p wave duration
       
 
data=p[1,:]
minimum = signal.argrelextrema(data, np.less)
 
if minimum[0][-1]>126 and minimum[0][-1]<135:               # trough the derivation
 
       pwd_query=minimum[0][-1]
       print(minimum)
 
else:
       startwert=p[1,0]                                    # through a threshold depending on the first value
      
       if startwert<-0.015:
          
            endwert=np.where((data>startwert+0.016))
            pwd_query=endwert[0][-1]
            
       elif startwert<-0.01:
          
             endwert=np.where((data>startwert+0.011))
             pwd_query=endwert[0][-1]
       
       elif startwert<0:
          
             endwert=np.where((data>startwert+0.004))
             pwd_query=endwert[0][-1]
       
       else:
             endwert=np.where((data>startwert+0.005))
             pwd_query=endwert[0][-1]
 
p_crop=p[:,0:pwd_query]
 
#find qt duration through a threshold
 
start=np.array(np.where(qrst_scaled[1,0:50]>0))
 
if start.size==0:
   
    start=np.array(np.where(qrst_scaled[1,0:50]>-0.01))
 
start_qrst=start[0][0]
end=np.array(np.where(qrst_scaled[1,200:-1]>0.0001))
qt_query= 200+end[0][-1]-start_qrst
 
'''
 
3. Resample Single Beat:
   
    -> resample the single beat ecg with new intervals for the poff_qon and rr intervals
   
'''
 
# generate a Multivariate Normal Distributions for the intervals (Poff-Qon, P wave duration, QT interval, RR interval)
mean_int = np.array([25.168874868559410,1.495199789695058e+02,4.019283911671924e+02,8.821697160883281e+02])
sigma_int = np.array([[237.652658986445,40.3315366745504,57.5015776851198,331.851675705798],
                     [40.3315366745504,162.246446266645,55.2128397137995,226.669795367321],
                     [57.5015776851198,55.2128397137995,1259.01117080808,3785.91814686185],
                     [331.851675705798,226.669795367321,3785.91814686185,24621.2574214757]])
 
# generate random samples from the intervals' multivariate distribution
 
rnd_int = np.random.multivariate_normal(mean_int, sigma_int, 10000000)
 
success=False
thres=3
while success==False:
   
    rnd_int = np.random.multivariate_normal(mean_int, sigma_int, 10000000)
    thres=thres+1
    rnd_qt=np.nonzero(abs(rnd_int[:,2]-qt_query)<thres)
    rnd_pwd= np.nonzero(abs(rnd_int[:,1]-pwd_query)<thres)
    val=np.intersect1d(rnd_qt,rnd_pwd)
      
    if val.size==0:
        success=False
    else:
        success=True
# take the poff_qon and rr intervals that correspondend to the duration of the qrst and p wave 
 
R_constraint=np.zeros([val.size,2])
R_constraint[:,0]=rnd_int[val,0]
R_constraint[:,1]=rnd_int[val,3]
 # fit a Gaussion mixture model to the data
GMModel_constraint_int=GaussianMixture(n_components=1, random_state=0).fit(R_constraint)
 
# the new poff_qon and rr intervals are the mean values of the Gaussian mixture model
 
poff_qon_int=round(GMModel_constraint_int.means_[0][0])
rr_int=round(GMModel_constraint_int.means_[0][1])
  
qt_orig=qrst_scaled[1,:].size
 
# resampled the single beat ECG with the newly calculated intervals:
 
ecg_beat=[0]*12
for leadID in range(0,12):
   
    m = (qrst_scaled[leadID,0]-p_crop[leadID,-1])/(poff_qon_int)             # slope for the newly calculated poff_qon interval
    pq_sig = p_crop[leadID,-1] + np.arange(1,poff_qon_int+1)*m               # straight line between end of the p wave and the beginning of the q wave over the newly calculated interval
  
    ecg_beat[leadID]=list(p_crop[leadID,:])
    ecg_beat[leadID].extend(pq_sig)                                          # connect the p wave and the straight line
    
    li4=qrst_scaled[leadID,:].tolist()
    ecg_beat[leadID].extend(li4)                                             # connect the qrst wave with the p wave and the straight line
       
    tp_dur = rr_int - (pwd_query+poff_qon_int+qt_orig)                       # calculate the interval between the end of the qrst complex and the bginning of the new ecg
    m = (p_crop[leadID,0]-qrst_scaled[leadID,-1])/(tp_dur)                   # slope for the newly calculated tp interval
    tp_sig = qrst_scaled[leadID,-1] + np.arange(1,tp_dur+1)*m                # straight line from the end of the t wave to the beginning of a new p wave
    li_tp=tp_sig.tolist()
    li5=tp_sig.tolist()
    
    ecg_beat[leadID].extend(tp_sig)
   
 
'''
 
4. Resample the time intervals between the single heart beats and the duration of the qrst complex
 
'''
# generate RR time series
 
from RRtimeserie import *
f_getRRtimeserie(rr_int*0.001,15)
tau= 1000*tau        
 
qt_mod=np.empty((15,1))*np.nan
tp_mod=np.empty((15,1))*np.nan
 # find qt and tp intervals for the RR time series
thres=3
test=False
 
for t_ID in range(0,15):
   
    rnd_int = np.random.multivariate_normal(mean_int, sigma_int, 10000000)
      
     #find all random numbers close to the simulated query values
      
    rnd_rr=np.nonzero(abs(rnd_int[:,3]-tau[t_ID,0])<thres)
    rnd_pqi=np.nonzero(abs(rnd_int[:,0]-poff_qon_int)<thres)
    rnd_pwd=np.nonzero(abs(rnd_int[:,1]-pwd_query)<thres)
    val=np.intersect1d(rnd_rr,rnd_pwd)
 
    while val.size<2:
           
            thres=thres+1
 
            rnd_rr=np.nonzero(abs(rnd_int[:,3]-tau[t_ID,0])<thres)
            rnd_pqi=np.nonzero(abs(rnd_int[:,0]-poff_qon_int)<thres)
            rnd_pwd=np.nonzero(abs(rnd_int[:,1]-pwd_query)<thres)
            val=np.intersect1d(rnd_rr,rnd_pwd)
    
    # take the qt intervals that correspond to the found rr intervals and p duration
   
    R_constraint=rnd_int[val,2]
    R_constraint=R_constraint.reshape(-1,1)
   
    # fit a Gaussion mixture model to the data
       
    GMModel_constraint_int=GaussianMixture(n_components=1).fit(R_constraint)
   
    # new qt interval for every tau
   
    qt_mod[t_ID,0]=round(GMModel_constraint_int.means_[0][0])
    qt_mod=qt_mod.astype(int)
   
    # new tp interval for every tau
      
    tp_mod[t_ID,0]=tau[t_ID,0]-poff_qon_int-pwd_query-qt_mod[t_ID,0]
 
# put everything togehter
# -> take the single beat and change QT segment as well as TP segment.
 
ecg_beats=[ [0] * 12  for t in range(15)]
 
for t_ID in range(0,15):
   
    for leadID in range (0,12):
       
        # resample the qrst segment with the new qt value
        rate=int(qt_mod[t_ID,0]/qt_query*len(qrst_scaled[leadID,:]))
        qt_new=signal.resample(qrst_scaled[leadID,:],rate)
       
        # new straight line for the tp interval
        tp_new=tau[t_ID,0]-poff_qon_int-pwd_query-len(qt_new)
        m = (p_crop[leadID,0]-qrst_scaled[leadID,-1])/(tp_new)
        tp_sig=qrst_scaled[leadID,-1]+np.arange(1,tp_new+1)*m
          
        li=ecg_beat[leadID][0:pwd_query+poff_qon_int-2]
        li1=qt_new.tolist()
        li2=tp_sig.tolist()
          
        li.extend(li1)
        li.extend(li2)
        
        ecg_beats[t_ID][leadID]=li
          
 # arrange the single heartbeats in a row            
            
e=list(map(list, zip(*ecg_beats)))
    
ecg_series=[ []  for t in range(12)]
for i in range(0,12):# iterating over the data
    for item in e[i]:
        # appending elements to the flat_list
        ecg_series[i] += item  
 
 
ecg_series_d=signal.decimate(ecg_series,2)
 
# plot unfiltered ECG
 
plt.plot(ecg_series[1][0:10000])
plt.title('Result:10s ECG')
plt.show()
 
# save the unfiltered ECG
 
name_ungefiltert=args.Speicherort+'/probe_ungefiltert.csv'
 
with open(name_ungefiltert, 'w') as file:
     writer = csv.writer(file, delimiter=',')
     writer.writerows(ecg_series_d[:][0:5000])
 
'''
5. Add noise and filter the signal
'''
 
# Add Noise:
    
data_noise=loadmat("DATA_noises_real.mat")
mixture_of_noise=data_noise['mixture_of_noises'][:][:]
start_noise=random.randrange(0,len(mixture_of_noise[1])-len(ecg_series[1]))
noise_amp = 0.03
 
noise=mixture_of_noise[0:12,start_noise:start_noise+len(ecg_series[1])]
 
ecg_series_noise=ecg_series+noise_amp*noise
 
ecg_noise_d=signal.decimate(ecg_series_noise,2)
 
 
# save the ECG with noise
 
name_ecg_noise=args.Speicherort+'/probe_with_noise.csv'
 
with open(name_ecg_noise, 'w') as file:
     writer = csv.writer(file, delimiter=',')
     writer.writerows(ecg_noise_d[:][0:5000])
 
 
# High and Lowpass Filer:
 
sos_h=signal.butter(3,0.5,btype='highpass',fs=1000,output='sos')
sos_l=signal.butter(3,150,btype='lowpass', fs=1000,output='sos')
 
 
filtered=signal.sosfiltfilt(sos_h,ecg_series_noise)
filtered_1=signal.sosfiltfilt(sos_l,filtered)
 
filtered_offset=[0]*11
for E in range(0,11):
 
    offset= -1 *filtered_1[E][0]
 
    filtered_offset[E]= filtered_1[E] + offset
 
 
ecg_d=signal.decimate(filtered_offset,2)
 
# save the filtered ECG
 
name=args.Speicherort+'/probe_gefiltert.csv'
 
with open(name, 'w') as file:
    writer = csv.writer(file, delimiter=',')
    writer.writerows(ecg_d)
 
