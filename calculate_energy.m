function Enrg=calculate_energy(d)
% Computes the energy per sample (average power) of the given detail coefficient
%
% INPUT:
% - d: detail coefficient 
% OUTPUT:
% - Enrg: energy per sample (average power)
%
% Original code by Erik Alonso

Enrg=1000*sum(d.^2)./length(d);