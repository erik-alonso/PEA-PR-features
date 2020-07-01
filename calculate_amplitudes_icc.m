function [mean_PT,var_PT,std_PT]=calculate_amplitudes_icc(icc,posp,fs,t_trans)
% Computes mean, variance and standard deviation of the peak-to-trough
% amplitude of the impedance circulation component
%
% INPUT:
% - icc: impedance circulation component
% - posp: instants in seconds of the QRS complexes detected by Hamilton-Tompkins (HT) detector 
% - fs: sampling rate of ICC
% - t_trans: initial signal interval (in seconds) not to be analyzed to avoid transient of the adaptive filter in the ICC 
% OUTPUT:
% - mean_PT: mean peak-to-trough amplitude of the ICC
% - var_PT : variance of the peak-to-trough amplitude of the ICC
% - std_PT : standard deviation of the peak-to-trough amplitude of the ICC
%
% Original code by Erik Alonso


cc_pp_array=[];
cc_max_array=[];
cc_min_array=[];
posp=posp(posp>t_trans);

for cc=1:length(posp)
    posp_actual=posp(cc);
    if cc==length(posp)
        fin=length(icc);
        ini=round(posp_actual.*fs);
    else
        posp_sig=posp(cc+1);
        ini=round(posp_actual.*fs);
        fin=round(posp_sig.*fs);
    end
    x=icc(ini:fin);
    max_p=[];
    max_v=[];
    for i=2:length(x)-1
        if (x(i)>x(i-1))&&(x(i)>x(i+1))
            max_p=[max_p i];
            max_v=[max_v x(i)];
        end
    end
    
    [cc_max_v,cc_max_p]=max(max_v);
    cc_max_p=max_p(cc_max_p);
    if isempty(max_v)
        [cc_max_v,cc_max_p]=max(x);
    end
    
    min_p=[];
    min_v=[];
    for i=2:length(x)-1
        if (x(i)<x(i-1))&&(x(i)<x(i+1))
            min_p=[min_p i];
            min_v=[min_v x(i)];
        end
    end
    
    [cc_min_v,cc_min_p]=min(min_v);
    cc_min_p=min_p(cc_min_p);
    if isempty(min_v)
        [cc_min_v,cc_min_p]=min(x);
    end
    
    if ~isempty(cc_max_v) && ~isempty(cc_min_v)
        cc_max_array=[cc_max_array ini+cc_max_p];
        cc_min_array=[cc_min_array ini+cc_min_p];
        cc_pp=cc_max_v-cc_min_v;
        cc_pp_array=[cc_pp_array cc_pp];
    end
end
mean_PT=mean(cc_pp_array);
var_PT=var(cc_pp_array);
std_PT=std(cc_pp_array);

    
