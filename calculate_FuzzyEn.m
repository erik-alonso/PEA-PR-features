function [FuzzyEn] = calculate_FuzzyEn(series,dim,r,n)
%Feature from:
% "Characterization of surface EMG signal based on fuzzy entropy"
% Chen et al.2007; IEEE Transactions on neural systems and rehabilitation
% engineering 15(2):266-272.
%
% "Measuring complexity using FuzzyEn, ApEn, and SampEn" 
% Chen et al.2009; Medical Engineering \& Physics 31(1):61-68.
%
% INPUT:
% - series: input time series
% - m: maximum template length
% - r:  matching tolerance[mV]
% - n : the step of the fuzzy exponential function
% OUTPUT:
% - FuzzyEn: Fuzzy Entropy.
%
% Original code by Jesús Monge Álvarez 
% Updated by U Irusta and B Chicote
%
% DATE: 21/02/2017
% Revised by Unai Irusta (unia.irusta@ehu.eus) and Beatriz Chicote (beatriz.chicote@ehu.eus)
% to reproduce the definition introduced by Chen et al. "Measuring complexity using FuzzyEn, ApEn, and SampEn" 
% Change is in the computation of similarity: 
%       simi = exp(((-1)*((dist).^n/r))); -> simi = exp(((-1)*((dist/r).^n)));

% Do not normalize time-series
% This is done according to Chicote et al 2016.....
% It is dependent on VF-amplitude
N = length(series);
phi = zeros(1,2);

%%% This is the orifinal code by Jesús Monge Álvarez
for j = 1:2
    m = dim+j-1; % 'm' is the embbeding dimension used each iteration
    % Pre-definition of the varialbes for computational efficiency:
    patterns = zeros(m,N-m+1);
    aux = zeros(1,N-m+1);
    
    % First, we compose the patterns
    % The columns of the matrix 'patterns' will be the (N-m+1) patterns of 'm' length:
    if m == 1 % If the embedding dimension is 1, each sample is a pattern
        patterns = series;
    else % Otherwise, we build the patterns of length 'm':
        for i = 1:m
            patterns(i,:) = series(i:N-m+i);
        end
    end
    % We substract the baseline of each pattern to itself:
    for i = 1:N-m+1
        patterns(:,i) = patterns(:,i) - (mean(patterns(:,i)));
    end

    % This loop goes over the columns of matrix 'patterns':
    for i = 1:N-m
        % Second, we compute the maximum absolut distance between the
        % scalar components of the current pattern and the rest:
        if m == 1 
            dist = abs(patterns - repmat(patterns(:,i),1,N-m+1));
        else
            dist = max(abs(patterns - repmat(patterns(:,i),1,N-m+1)));
        end
       % Third, we get the degree of similarity:
       %%% Changed JMA code to the modified definition used by Chen et al 2009.
       simi = exp(((-1)*((dist/r).^n)));
       % We average all the degrees of similarity for the current pattern:
       aux(i) = (sum(simi)-1)/(N-m-1); % We substract 1 to the sum to avoid the self-comparison
    end

    % Finally, we get the 'phy' parameter as the as the mean of the first
    % 'N-m' averaged drgees of similarity:
    phi(j) = sum(aux)/(N-m);
end
FuzzyEn = log(phi(1)) - log(phi(2));
