function MA=calculate_area_icc(icc,t_trans,fs)
% Computes the mean area of the ICC
%
% INPUT:
% - icc: impedance circulation component 
% - t_trans: initial signal interval (in seconds) not to be analyzed to avoid transient of the adaptive filter in the ICC 
% - fs: sampling rate of ICC
% OUTPUT:
% - MA: Mean Area of the ICC
%
% Original code by Erik Alonso

cc=icc(t_trans*fs+1:end);
MA=sum(abs(cc))/length(cc);