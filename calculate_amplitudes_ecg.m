function [mean_PPA, var_PPA, std_PPA]=calculate_amplitudes_ecg(ecg,posp,fs)
% Computes mean, variance and standard deviation of the peak-to-peak
% amplitude of the QRS complexes
%
% INPUT:
% - ecg: preprocessed/denoised ECG
% - posp: instants in seconds of the QRS complexes detected by Hamilton-Tompkins (HT) detector 
% - fs: sampling rate of ECG
% OUTPUT:
% - mean_PPA: mean peak-to-peak amplitude of the QRS complexes
% - var_PPA:  variance of the peak-to-peak amplitude of the QRS complexes
% - std_PPA:  standard deviation of the peak-to-peak amplitude of the QRS complexes
%
% Original code by Erik Alonso

ecg_pp_array=[];
ecg_max_array=[];
ecg_min_array=[];
posp_int=posp(2:end)-posp(1:end-1);
interval=round(mean(posp_int)/2*fs);


for cc=1:length(posp)
    posp_actual=posp(cc);
    if length(posp)>1
        if cc==1
            interval2=round((posp(cc+1)-posp_actual)/2*fs);
            ini=max(1,posp_actual*fs-interval);
            fin=posp_actual*fs+interval2;
        elseif cc==length(posp)
            interval1=round((posp_actual-posp(cc-1))/2*fs);
            fin=min(length(ecg),posp_actual*fs+interval);
            ini=posp_actual*fs-interval1;
        else
            interval1=round((posp_actual-posp(cc-1))/2*fs);
            interval2=round((posp(cc+1)-posp_actual)/2*fs);
            ini=posp_actual*fs-interval1;
            fin=posp_actual*fs+interval2;
        end
    else %Just a single QRS complex
        th=round(fs/4);
        ini=max(1,posp_actual*fs-th);
        fin=min(length(ecg),posp_actual*fs+th);
    end
    ini=int64(ini);
    fin=int64(fin);
    x=ecg(ini:fin);
    max_p=[];
    max_v=[];
    for i=2:length(x)-1
        if (x(i)>x(i-1))&&(x(i)>x(i+1))
            max_p=[max_p i];
            max_v=[max_v x(i)];
        end   
    end
        
    [ecg_max_v,ecg_max_p]=max(max_v);
    ecg_max_p=max_p(ecg_max_p);
    
    min_p=[];
    min_v=[];
    for i=2:length(x)-1
        if (x(i)<x(i-1))&&(x(i)<x(i+1))
            min_p=[min_p i];
            min_v=[min_v x(i)];
        end   
    end

   [ecg_min_v,ecg_min_p]=min(min_v);
   ecg_min_p=min_p(ecg_min_p);
   if ~isempty(ecg_max_v) && ~isempty(ecg_min_v)
    ecg_max_array=[ecg_max_array ini+ecg_max_p];
    ecg_min_array=[ecg_min_array ini+ecg_min_p];
    ecg_pp=ecg_max_v-ecg_min_v;
    ecg_pp_array=[ecg_pp_array ecg_pp];
   end
end

mean_PPA=mean(ecg_pp_array);
var_PPA=var(ecg_pp_array);
std_PPA=std(ecg_pp_array);
