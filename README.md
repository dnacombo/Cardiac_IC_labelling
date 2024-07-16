# Cardiac_IC_labelling
Find automatically cardiac IC (Independent Component) without ECG

% Author
% % % % % % % % % % % % % % % % % 
% Pierre Champetier (2024)
% Contact: pi.champetier@gmail.com

% Aim
% % % % % % % % % % % % % % % % %
% This function takes EEG independent components (IC) as inputs, and finds the cardiac IC. 
% --> To do so:
% 1) The script detects cardiac events (PQRST) in all IC.
% 2) It computes several features (std/mean ratio, skewness, kurtosis...) for the distribution of RR intervals, R, Q and S amplitudes. 
% 3) It identifies the top 3 IC with the lowest values for each feature (or highest values depending on the feature), and gives a score (1 or 0) for the feature.
% 4) It computes a global score (sum of score/nb of feature).
% 5) After sanity check (physiological beat per minute (bpm) and regularity of IC timecourse), cardiac IC are identified based on the global score:
%   -method 1: cardiac IC if global score >= threshold (ex: 0.5)
%   -method 2: cardiac if global score >= mean(all_global scores) + n*std(all_global_scores)


% Inputs
% % % % % % % % % % % % % % % % % 
% 1) comp --> your IC in FieldTrip format
%     Ex: EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',round(dataRank/5)); % Run your ICA with EEGLAB
%         comp = eeglab2fieldtrip(EEG, 'comp'); % Convert into FieldTrip format
%
% 2) method_chosen --> 'absolute_threshold' (method 1) or 'mean_std' (method 2)
% 3) plot_heart_IC --> 1 or 0 (to plot the IC labelled as cardiac)
% 4) path_output --> path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
% 5) file_info --> name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)


% Usage
% % % % % % % % % % % % % % % % % 
% Minimum inputs:
% If plot_heart_IC == 0
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(comp, method_chosen, plot_heart_IC); 
%
% If plot_heart_IC == 1, you can also give your output path and the name of your recording to have it in the file name of your plot
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(comp, method_chosen, plot_heart_IC, path_output, file_info);


% Parameters
% % % % % % % % % % % % % % % % % 
% You can easily change the parameters used at the begining of the function "A_fct_find_cardiac_IC.m"
% nb_IC_wanted = 3; % The heart IC must be in the top 3 of all IC for skewness, kurtosis... of Rampl, RRintervals...
% bpm_min = 45;
% bpm_max = 90;
% threshold_cond_IC_method1 = 0.5;
% threshold_std_method2 = 2.5;
% threshold_regularity_signal_minmax = 1.5;
% min_recording_duration_sec = 20;


% Toolbox used to detect cardiac events
% % % % % % % % % % % % % % % % % 
% Based on R. Sanghavi, F. Chheda, S. Kanchan and S. Kadge, "Detection Of Atrial Fibrillation in Electrocardiogram Signals using Machine Learning," 2021 2nd Global Conference for Advancement in Technology (GCAT), 2021, pp. 1-6, doi: 10.1109/GCAT52182.2021.9587664.
% Info : https://fr.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox
% Citation pour cette source: Rohan Sanghavi (2024). ECG SIGNAL PQRST PEAK DETECTION TOOLBOX (https://www.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox), MATLAB Central File


% Dependencies
% % % % % % % % % % % % % % % % % % 
% 1) A_fct_test_unif.m (To test the uniform distribution of cardiac events)
% 2) ECG_PQRST_VERSION_3 (Toolbox used to detect cardiac events) 
% --> Rq: I did small modifications of 'compute_fudicial_peaks_live17_c' and 'preprocess_window_ecg' functions to avoid errors with some IC (all modifications are preceeded by a comment "PIERRE").

