function [mean_RR, std_RR, nQRS]=calculate_RR_features(posp)
% Computes mean and standard deviation of the RR interval and number of QRS
% complexes in the analyzed segment
%
% INPUT:
% - posp: instants in seconds of the QRS complexes detected by Hamilton-Tompkins (HT) detector 
% OUTPUT:
% - mean_PPA: mean peak-to-peak amplitude of the QRS complexes
% - var_PPA:  variance of the peak-to-peak amplitude of the QRS complexes
% - std_PPA:  standard deviation of the peak-to-peak amplitude of the QRS complexes
%
% Original code by Erik Alonso

posp_int=posp(2:end)-posp(1:end-1);
mean_RR=mean(posp_int);
std_RR=std(posp_int);
nQRS=length(posp);