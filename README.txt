In this folder the user can find the MATLAB scripts used in the following work:
Alonso et al. A Machine Learning Framework for Pulse Detection During Out-of-hospital
Cardiac Arrest
to compute the features to discriminate between PEA and PR rhythms. Briefly:

1) calculate_amplitudes_ecg.m
Script to compute mean, variance and standard deviation of the peak-to-peak
amplitude of the QRS complexes.

2) calculate_amplitudes_icc.m
Script to compute the mean, variance and standard deviation of the peak-to-trough
amplitude of the impedance circulation component.

3) calculate_area_icc.m
Script to compute the mean area of a signal.

4) calculate_crossPower.m
Script to compute the cross-power between two signals.

5) calculate_energy.m
Script to compute the energy per sample (average power) of the given signal.

6) calculate_FuzzyEn.m
Script to compute the Fuzzy Entropy of a given signal.

7) calculate_Kurtosis.m
Script to compute the Kurtosis of a given signal.

8) calculate_QRS_width.m
Script to compute the median, mean and standard deviation of QRS complex width

9) calculate_RR_features.m
Script to compute the mean and standard deviation of the RR interval and number of QRS
complexes in the analyzed segment

10) calculate_spectral_measures.m
Script to compute AMSA (AMplitude Spectrum Area) of preprocessed/denoised ECG and its energy 
in the band 17.5 - 30 Hz.
