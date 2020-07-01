function [AMSA,Sxx]=calculate_spectral_measures(ecg,fs,tw)
% Computes AMSA (AMplitude Spectrum Area) of preprocessed/denoised ECG as in
% Ristagno et al. Amplitude Spectrum Area to Guide Defibrillation.
% Validation on 1617 Patients With Ventricular Fibrillation.
% Circulation 2015; 131: 478-487.
%
% Computes also energy in the band 17.5 - 30 Hz as in 
% Jekova et al. Real time detection of ventricular fibrillation and
% tachycardia. Physiol. Meas. 2004; 25: 1167–1178.
%
% INPUT:
% - ecg: preprocessed/denoised ECG
% - fs: sampling rate of ECG
% - tw: duration of analysis window (in seconds)
% OUTPUT:
% - AMSA: AMplitude Spectrum Area
% - Sxx:  Energy in the 17.5 - 30 Hz frequency band
%
% Original code by Erik Alonso

ecg=ecg(1:round(tw*fs));
w=tukeywin(length(ecg),0.2)';
ecg_w=ecg.*w;

n=nextpow2(length(ecg_w));
NFFT=max(2^n,4096);

[XW,F]=freqz(ecg_w,1,NFFT,fs);

% AMSA
XW=abs(XW);f1=2;f2=30;
AMSA=sum(XW(F>f1 & F<f2).*F(F>f1 & F<f2))/NFFT;

% Energy in the 17.5 - 30 Hz frequency band
Pxx=XW.^2;fc1=17.5;fc2=30;
Pxx=sum(Pxx(F>fc1 & F<fc2));
Sxx=Pxx/NFFT*fs/2;