# Cardiac_IC_labelling
Find automatically cardiac IC (Independent Component) without ECG


%%%%%% AUTHOR %%%%%%

Pierre Champetier (2024)
Contact: pi.champetier@gmail.com

%%%%%% MAIN FUNCTION %%%%%%

A_fct_find_cardiac_IC.m


%%%%%% AIM %%%%%%

This function takes EEG independent components (IC) as inputs, and finds the cardiac IC. 
  --> To do so:
1) The script detects cardiac events (PQRST) in all IC.
2) It computes several features (std/mean ratio, skewness, kurtosis...) for the distribution of RR intervals, R, Q and S amplitudes. 
3) It identifies the top 3 IC with the lowest values for each feature (or highest values depending on the feature), and gives a score (1 or 0) for the feature.
4) It computes a global score (sum of score/nb of feature).
5) After sanity check (physiological beat per minute (bpm) and regularity of IC timecourse), cardiac IC are identified based on the global score:
  -method 1: cardiac IC if global score >= threshold (ex: 0.5)
  -method 2: cardiac if global score >= mean(all_global scores) + n*std(all_global_scores)



%%%%%% INPUTS %%%%%%

1) cfg --> all the parameters including:
     -method_chosen --> 'absolute_threshold' (method 1) or 'mean_std' (method 2)
     -plot_heart_IC --> 1 or 0 (to plot the IC labelled as cardiac)
     -path_output --> path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
     -file_info --> name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
     -nb_IC_wanted --> number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
     -bpm_min and bpm_max  --> expected heart beat per min, for sanity check [default: 45 and 90]
     -threshold_cond_IC_method1 --> minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
     -threshold_std_method2 --> if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
     -min_recording_duration_sec --> minimum duration (in sec) of the IC timecourse (default: 20]
     -mini_bouts_duration_for_SignalAmplRange --> for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
     -threshold_regularity_signal_minmax --> For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]

2) comp --> your IC in FieldTrip format 



%%%%%% USAGE%%%%%%

[rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

Example 1 (all default parameters):
EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',round(dataRank/5)); %Run your ICA with EEGLAB
comp = eeglab2fieldtrip(EEG, 'comp'); % Convert into FieldTrip format
cfg = [];
[rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

Example 2 (If you want to use other parameters, specify them in cfg):
cfg = [];
cfg.nb_IC_wanted = 5; % The heart IC must be in the top 5 of all IC for skewness, kurtosis... of Rampl, RRintervals...
cfg.bpm_max = 120; % Max physiological bpm
[rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

Rq: If plot_heart_IC == 1, you must specify cfg.path_output and cfg.file_info (output path and the name of your recording to have it in the file name of your plot) 



%%%%%% TOOLBOX USED TO DETECT CARDIAC EVENTS %%%%%%

Based on R. Sanghavi, F. Chheda, S. Kanchan and S. Kadge, "Detection Of Atrial Fibrillation in Electrocardiogram Signals using Machine Learning," 2021 2nd Global Conference for Advancement in Technology (GCAT), 2021, pp. 1-6, doi: 10.1109/GCAT52182.2021.9587664.
Info : https://fr.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox
Citation pour cette source: Rohan Sanghavi (2024). ECG SIGNAL PQRST PEAK DETECTION TOOLBOX (https://www.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox), MATLAB Central File



%%%%%% DEPENDENCIES %%%%%%

1) A_fct_test_unif.m (To test the uniform distribution of cardiac events)
2) ECG_PQRST_VERSION_3 (Toolbox used to detect cardiac events) 
--> Rq: I did small modifications of 'compute_fudicial_peaks_live17_c' and 'preprocess_window_ecg' functions to avoid errors with some IC (all modifications are preceeded by a comment "PIERRE").

