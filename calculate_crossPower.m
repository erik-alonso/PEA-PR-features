function CP=calculate_crossPower(ecg,icc,fs_ecg,fs_icc,t_trans)
%Computes the cross-power between preprocessed ECG and ICC as in
%Ruiz et al. Circulation assessment by automated external defibrillators during
%cardiopulmonary resuscitation. Resuscitation 2018; 128: 158-163.
%
% INPUT:
% - ecg: preprocessed/denoised ECG
% - icc: impedance circulation component
% - fs_ecg: sampling rate of ECG
% - fs_icc: sampling rate of ICC
% - t_trans: initial signal interval (in seconds) not to be analyzed to avoid transient of the adaptive filter in the ICC 
% OUTPUT:
% - CP: Cross Power
%
% Original code by Erik Alonso

icc=icc.*1000;%From Ohm to mOhm
% Evaluate if resample is needed
if fs_ecg>fs_icc
    nmf=100;
    ecg1=[ones(1,100)*ecg(1) ecg ones(1,100)*ecg(end)];
    ecg2=resample(ecg1,fs_icc,fs_ecg);
    ecg=ecg2(nmf*fs_icc/fs_ecg+1:end-nmf*fs_icc/fs_ecg);
elseif fs_ecg<fs_icc
    nmf=100;
    icc1=[ones(1,100)*icc(1) icc ones(1,100)*icc(end)];
    icc2=resample(icc1,fs_ecg,fs_icc);
    icc=icc2(nmf*fs_ecg/fs_icc+1:end-nmf*fs_ecg/fs_icc);
end
ini=round(t_trans*fs_icc);% To avoid C of the adaptive extraction of the ICC
ecg=ecg(ini+1:end);
icc=icc(ini+1:end);
L=min(length(ecg),length(icc));
N=floor(L/2);
if N<L
    ecgw=ecg(1:N);
    iccw=icc(1:N);
    
    % Divide into two subsegments
    aCP(1)=2/N*sum(abs(ecgw(1:round(N/2))).*abs(iccw(1:round(N/2))));
    aCP(2)=2/N*sum(abs(ecgw(round(N/2)+1:end)).*abs(iccw(round(N/2)+1:end)));
    
else
    disp('Analysis window longer than segment length');
end

CP=min(aCP);