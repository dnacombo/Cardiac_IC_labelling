function [heart_IC, table_cardiac_IC, aaa_parameters_find_heart_IC] = fct_find_cardiac_IC(cfg, comp)

%% DOCSTRING

%%%%%% AUTHOR %%%%%%
% Pierre Champetier (2024) Contact: pi.champetier@gmail.com


% %%%%%% AIM %%%%%%
% This function takes EEG independent components (IC) as inputs, and finds the cardiac IC. --> To do so:
% 
% The script detects cardiac events (PQRST) in all IC.
% It computes several features (std/mean ratio, skewness, kurtosis...) for the distribution of RR intervals, R, Q and S amplitudes.
% It identifies the top 3 IC with the lowest values for each feature (or highest values depending on the feature), and gives a score (1 or 0) for the feature.
% It computes a global score (sum of score/nb of feature).
% After sanity check (physiological beat per minute (bpm) and regularity of IC timecourse), cardiac IC are identified based on the global score: -method 1: cardiac IC if global score >= threshold (ex: 0.5) -method 2: cardiac if global score >= mean(all_global scores) + n*std(all_global_scores)


% %%%%%% INPUTS %%%%%%
% 1) cfg --> all the parameters including:
    % -method_chosen --> 'absolute_threshold' (method 1) or 'mean_std' (method 2)
    % -plot_heart_IC --> 1 or 0 (to plot the IC labelled as cardiac)
    % -path_output --> path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
    % -file_info --> name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
    % -nb_IC_wanted --> number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
    % -bpm_min and bpm_max  --> expected heart beat per min, for sanity check [default: 45 and 90]
    % -threshold_cond_IC_method1 --> minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
    % -threshold_std_method2 --> if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
    % -min_recording_duration_sec --> minimum duration (in sec) of the IC timecourse (default: 20]
    % -mini_bouts_duration_for_SignalAmplRange --> for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
    % -threshold_regularity_signal_minmax --> For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]

% 2) comp --> your IC in FieldTrip format 


% %%%%%% USAGE %%%%%%
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);
%
% Example 1 (all default parameters):
% EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',round(dataRank/5)); %Run your ICA with EEGLAB
% comp = eeglab2fieldtrip(EEG, 'comp'); % Convert into FieldTrip format
% cfg = [];
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

% Example 2 (If you want to use other parameters, specify them in cfg):
% cfg = [];
% cfg.nb_IC_wanted = 5; % The heart IC must be in the top 5 of all IC for skewness, kurtosis... of Rampl, RRintervals...
% cfg.bpm_max = 120; % Max physiological bpm
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

% Rq: If plot_heart_IC == 1, you must specify cfg.path_output and cfg.file_info (output path and the name of your recording to have it in the file name of your plot) 


% %%%%%% TOOLBOX USED TO DETECT CARDIAC EVENTS %%%%%%
% Based on R. Sanghavi, F. Chheda, S. Kanchan and S. Kadge, "Detection Of Atrial Fibrillation in Electrocardiogram Signals using Machine Learning," 2021 
% 2nd Global Conference for Advancement in Technology (GCAT), 2021, pp. 1-6, doi: 10.1109/GCAT52182.2021.9587664. 
% Info : https://fr.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox 
% Citation pour cette source: Rohan Sanghavi (2024). ECG SIGNAL PQRST PEAK DETECTION TOOLBOX 
% (https://www.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox), MATLAB Central File


% %%%%%% DEPENDENCIES %%%%%%
% A_fct_test_unif.m (To test the uniform distribution of cardiac events)
% ECG_PQRST_VERSION_3_PC (Toolbox used to detect cardiac events) --> Rq: I did small modifications of 'compute_fudicial_peaks_live17_c' and 'preprocess_window_ecg' functions to avoid errors with some IC (all modifications are preceeded by a comment "PIERRE").



%% Extract parameters

% method_chosen
if isfield(cfg, 'method_chosen') == 1
    method_chosen = cfg.method_chosen;
else
    method_chosen = 'absolute_threshold'; % absolute_threshold (method 1) or mean_std (method 2)
end

% plot_heart_IC
if isfield(cfg, 'plot_heart_IC') == 1
    plot_heart_IC = cfg.plot_heart_IC;
else
    plot_heart_IC = 0;
end

% path_output and file_info
if plot_heart_IC == 1
    try
        path_output = cfg.path_output;
        file_info = cfg.file_info;
    catch
        error('When plot_heart_IC == 1, you must specified cfg.path_output and cfg.file_info')
    end
end

% nb_IC_wanted
if isfield(cfg, 'nb_IC_wanted') == 1
    nb_IC_wanted = cfg.nb_IC_wanted;
else
    nb_IC_wanted = 3; % default
end

% bpm_min
if isfield(cfg, 'bpm_min') == 1
    bpm_min = cfg.bpm_min;
else
    bpm_min = 45;
end

% bpm_max
if isfield(cfg, 'bpm_max') == 1
    bpm_max = cfg.bpm_max;
else
    bpm_max = 90;
end

% threshold_cond_IC_method1
if isfield(cfg, 'threshold_cond_IC_method1') == 1
    threshold_cond_IC_method1 = cfg.threshold_cond_IC_method1;
else
    threshold_cond_IC_method1 = 0.5;
end

% threshold_std_method2
if isfield(cfg, 'threshold_std_method2') == 1
    threshold_std_method2 = cfg.threshold_std_method2;
else
    threshold_std_method2 = 2.5;
end

% threshold_regularity_signal_minmax
if isfield(cfg, 'threshold_regularity_signal_minmax') == 1
    threshold_regularity_signal_minmax = cfg.threshold_regularity_signal_minmax;
else
    threshold_regularity_signal_minmax = 1.5;
end

% min_recording_duration_sec
if isfield(cfg, 'min_recording_duration_sec') == 1
    min_recording_duration_sec = cfg.min_recording_duration_sec;
else
    min_recording_duration_sec = 20;
end

% mini_bouts_duration_for_SignalAmplRange
if isfield(cfg, 'mini_bouts_duration_for_SignalAmplRange') == 1
    mini_bouts_duration_for_SignalAmplRange = cfg.mini_bouts_duration_for_SignalAmplRange;
else
    mini_bouts_duration_for_SignalAmplRange = 10;
end


%% Extract parameters

aaa_parameters_find_heart_IC = [];
aaa_parameters_find_heart_IC.nb_IC_wanted = nb_IC_wanted;
aaa_parameters_find_heart_IC.bpm_min = bpm_min;
aaa_parameters_find_heart_IC.bpm_max = bpm_max;
aaa_parameters_find_heart_IC.threshold_regularity_signal_minmax = threshold_regularity_signal_minmax;
aaa_parameters_find_heart_IC.threshold_cond_IC_method1 = threshold_cond_IC_method1;
aaa_parameters_find_heart_IC.threshold_std_method2 = threshold_std_method2;
aaa_parameters_find_heart_IC.min_recording_duration_sec = min_recording_duration_sec;
aaa_parameters_find_heart_IC.mini_bouts_duration_for_SignalAmplRange = mini_bouts_duration_for_SignalAmplRange;

fs = comp.fsample;

% Checking that recording_duration_sec > mini_bouts_duration_for_SignalAmplRange
sample_tot = 0;
for i = 1:length(comp.trial)
    sample_tot = sample_tot + length(comp.trial{1,i});
end
recording_duration_sec = sample_tot / fs;

if recording_duration_sec < mini_bouts_duration_for_SignalAmplRange
    error('Error: attempt to divide the recording into bouts with longer than the original recording (for the Signal Ampl Range checking). Change the value of mini_bouts_duration_for_SignalAmplRange variable in the function.')
end


%% Extract timecourse of each IC for all 1s-mini trials

IC_timecourse = [];
for i = 1:size(comp.trial,2)
    IC_timecourse = [IC_timecourse, comp.trial{1,i}];
end


%% Plot IC
% for comp_iter = 1:5
% figure
% plot([1:length(IC_timecourse)], IC_timecourse(comp_iter,:))
% title(num2str(comp_iter))
% end


%% DETECT CARDIAC EVENTS

plot_data = 0; % (REFER FUNCTION DOCUMENTATION)

IC_bug = [];

% RR features
RR_mean_all = [];
RR_std_all = [];
RR_skew_all = [];
RR_kurt_all = [];
RR_ratio_std_mean_all_IC = [];
RR_chi2stat_all_IC = [];

bpm_all_IC = [];


% Rampl features
Rampl_kurt_all_IC = [];
Rampl_skew_all_IC = [];
Rampl_ratio_std_mean_all_IC = [];
Rampl_chi2stat_all_IC = [];

SignalAmpl_range_all_IC = [];

% Qampl features
Qampl_kurt_all_IC = [];
Qampl_skew_all_IC = [];
Qampl_ratio_std_mean_all_IC = [];
Qampl_chi2stat_all_IC = [];

Sampl_kurt_all_IC = [];
Sampl_skew_all_IC = [];
Sampl_ratio_std_mean_all_IC = [];
Sampl_chi2stat_all_IC = [];


for comp_iter = 1:length(comp.label)
    % Select IC iter
    ecg_1 = IC_timecourse(comp_iter,:);

    % Detect cardiac events
    try
        ecg_f = preprocess_window_ecg(ecg_1, fs);
        [locs_P,locs_Q,locs_R,locs_S,locs_T] = compute_fudicial_peaks_live17_c(ecg_f, fs, plot_data); % Why the nb of samples in ecg_f is lower than in ecg_1?? For some IC it's super low....

    catch
        IC_bug = [IC_bug, comp_iter]; % I don't know why but the function preprocess_window_ecg bug for some recordings
        continue;
    end


    % Check enough data
    duration_sec_data_iter = length(ecg_f) / fs;
    if duration_sec_data_iter < min_recording_duration_sec
        IC_bug = [IC_bug, comp_iter];
        continue
    end


    %% Plot cardiac events
    %
    % figure(1);
    % hold on
    % plot(ecg_f,'linewidth',1.5);
    % plot(locs_R,ecg_f(locs_R),'rv','MarkerFaceColor','r');
    % plot(locs_Q,ecg_f(locs_Q),'rs','MarkerFaceColor','g');
    % plot(locs_S,ecg_f(locs_S),'o','MarkerFaceColor','y');
    % plot(locs_P,ecg_f(locs_P),'v','MarkerFaceColor','m');
    % plot(locs_T,ecg_f(locs_T),'s','MarkerFaceColor','k');
    % % plot(T2,'linewidth',2,'Linestyle','--','Color','red');
    % % plot(T1,'linewidth',2,'Linestyle','--','Color','green');
    % legend('ECG Signal','R peaks','Q Peaks','S Peaks','P peaks','T peaks','T-P1','T-P2');
    % hold off
    % title('QRS complex being marked on the fitered ECG signal along with the P and T peaks');
    % xlabel('Samples');
    % ylabel('Normalised Amplitude');
    % grid on

    %% Characteristics of RR intervals

    % 1) std/mean ratio, skew, kurt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Computation
    % % % % % % %

    % All parameters available --> could be also used to improve the nb of conditions to check and then improve performance?
    % [mean_QR, mean_SS, mean_ST, mean_RS, std_QR, std_SS, std_ST, skew_QR, skew_SS, skew_ST, skew_RS, kurt_QR, kurt_SS, kurt_ST, kurt_RS, std_RS] = compute_mean_intervals2_c(locs_R,locs_S,locs_Q,locs_T,fs);
    % [RR_mean, QRS_mean, mean_QT, mean_PR, RR_std, RR_skew, RR_kurt, QRS_std, QRS_skew, QRS_kurt, std_PR, skew_PR, kurt_PR, std_QT, skew_QT, kurt_QT]= compute_mean_interval_c(locs_P,locs_Q,locs_R,locs_S,locs_T,fs);

    % Those we need
    [RR_mean, RR_std, RR_skew, RR_kurt]= compute_mean_interval_c(locs_P,locs_Q,locs_R,locs_S,locs_T,fs);
    bpm = 60 * length(locs_R)/(duration_sec_data_iter);


    % Gather with other IC
    % % % % % % % % % % % %
    RR_mean_all = [RR_mean_all; comp_iter, RR_mean];
    RR_std_all = [RR_std_all; comp_iter, RR_std];
    RR_ratio_std_mean_all_IC = [RR_ratio_std_mean_all_IC; comp_iter, RR_std/RR_mean];
    RR_skew_all = [RR_skew_all; comp_iter, RR_skew];
    RR_kurt_all = [RR_kurt_all; comp_iter, RR_kurt];
    bpm_all_IC = [bpm_all_IC; comp_iter, bpm];



    % 2) Test uniform distribution of RR intervals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(locs_R) == 0
        chi2stat = A_fct_test_unif(locs_R, fs);
    else
        chi2stat = NaN;
    end
    RR_chi2stat_all_IC = [RR_chi2stat_all_IC; comp_iter, chi2stat];



    %% Characteristics of the AMPLITUDE of R, S and Q peaks

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                            R PEAKS
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % 1) std/mean ratio, skew and kurt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rampl_ratio_std_mean_all_IC = [Rampl_ratio_std_mean_all_IC; comp_iter, std(ecg_f(locs_R)) / mean(ecg_f(locs_R))];
    Rampl_kurt_all_IC = [Rampl_kurt_all_IC; comp_iter, kurtosis(ecg_f(locs_R))];
    Rampl_skew_all_IC = [Rampl_skew_all_IC; comp_iter, skewness(ecg_f(locs_R))];

    % 2) Test uniform distribution of amplitude of R peaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(locs_R) == 0
        chi2stat = A_fct_test_unif(ecg_f(locs_R), fs);
    else
        chi2stat = NaN;
    end
    Rampl_chi2stat_all_IC = [Rampl_chi2stat_all_IC; comp_iter, chi2stat];




    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                            Q PEAKS
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % 1) std/mean ratio, skew and kurt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qampl_ratio_std_mean_all_IC = [Qampl_ratio_std_mean_all_IC; comp_iter, std(ecg_f(locs_Q)) / mean(ecg_f(locs_Q))];
    Qampl_kurt_all_IC = [Qampl_kurt_all_IC; comp_iter, kurtosis(ecg_f(locs_Q))];
    Qampl_skew_all_IC = [Qampl_skew_all_IC; comp_iter, skewness(ecg_f(locs_Q))];

    % 2) Test uniform distribution of amplitude of R peaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(locs_R) == 0   
        chi2stat = A_fct_test_unif(ecg_f(locs_Q), fs);
    else
        chi2stat = NaN;
    end
    Qampl_chi2stat_all_IC = [Qampl_chi2stat_all_IC; comp_iter, chi2stat];



    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                            S PEAKS
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % 1) std/mean ratio, skew and kurt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sampl_ratio_std_mean_all_IC = [Sampl_ratio_std_mean_all_IC; comp_iter, std(ecg_f(locs_S)) / mean(ecg_f(locs_S))];
    Sampl_kurt_all_IC = [Sampl_kurt_all_IC; comp_iter, kurtosis(ecg_f(locs_S))];
    Sampl_skew_all_IC = [Sampl_skew_all_IC; comp_iter, skewness(ecg_f(locs_S))];

    % 2) Test uniform distribution of amplitude of R peaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(locs_R) == 0   
        chi2stat = A_fct_test_unif(ecg_f(locs_S), fs);
    else
        chi2stat = NaN;
    end
    Sampl_chi2stat_all_IC = [Sampl_chi2stat_all_IC; comp_iter, chi2stat];



    %% Metric to check homogeneous SignalAmpl across the recording
    % Divide the recordings into mini segments and check that the abs(mean)
    % is similar for each segment

    nbr_mini_segments = round((length(IC_timecourse(comp_iter,:)) / fs) / mini_bouts_duration_for_SignalAmplRange);

    % Extract the absolute value of the signal
    array = abs(IC_timecourse(comp_iter,:));

    % Determine the size of each small array
    sizeOfSmallArray = floor(length(array)/ nbr_mini_segments);

    % Preallocate cell array for the small arrays
    smallArrays = cell(1, nbr_mini_segments);

    % Fill the small arrays
    for i = 1:nbr_mini_segments
        startIdx = (i-1) * sizeOfSmallArray + 1;
        if i < nbr_mini_segments
            endIdx = i * sizeOfSmallArray;
        else
            % Last array takes any remaining elements
            endIdx = length(array);
        end
        smallArrays{i} = array(startIdx:endIdx);
    end

    % Compute mean for each segment
    SignalAmpl_mean_segments = [];
    for i = 1:nbr_mini_segments
        SignalAmpl_mean_segments = [SignalAmpl_mean_segments, mean(smallArrays{i})];
    end

    % Compute a metric to evaluate if the signal is regular over the recording (max-min)/mean
    SignalAmpl_range = (max(SignalAmpl_mean_segments) - min(SignalAmpl_mean_segments)) / min(SignalAmpl_mean_segments);
    SignalAmpl_range_all_IC = [SignalAmpl_range_all_IC; comp_iter, SignalAmpl_range];

end




%% Find heart IC

% 1) Manual visualization
%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(Rampl_kurt_all_IC(:,1), Rampl_kurt_all_IC(:,2))
% figure
% plot(RR_ratio_std_mean_all_IC(:,1), RR_ratio_std_mean_all_IC(:,2))
% figure
% plot(Rampl_skew_all_IC(:,1), Rampl_skew_all_IC(:,2))

% [sortedAbsArray, sortedIndices] = sort(abs(skew_RR_all_IC(:,2)));


% 2) Analyze distribution of RR_intervals and amplitude of R, S and Q peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the percentile we need to select the nb of IC chosen by user
nbr_IC_total = length(comp.label) - length(IC_bug);
perc = round(nb_IC_wanted * 100 /nbr_IC_total);


% RR_intervals
idx_score_std_mean = find(RR_ratio_std_mean_all_IC(:,2) < prctile(RR_ratio_std_mean_all_IC(:,2), perc)  );
IC_RR_std_mean = RR_ratio_std_mean_all_IC(idx_score_std_mean,1);

idx_score_skew = find(RR_skew_all(:,2) < prctile(RR_skew_all(:,2), perc)  );
IC_RR_skew = RR_skew_all(idx_score_skew,1);

idx_score_kurt = find(RR_kurt_all(:,2) < prctile(RR_kurt_all(:,2), perc)  );
IC_RR_kurt = RR_kurt_all(idx_score_kurt,1);

idx_chi2stat = find(RR_chi2stat_all_IC(:,2) < prctile(RR_chi2stat_all_IC(:,2), perc)  );
IC_RR_chi2stat = RR_chi2stat_all_IC(idx_chi2stat,1);


% Amplitude of R peaks
idx_score_RR = find(Rampl_ratio_std_mean_all_IC(:,2) < prctile(Rampl_ratio_std_mean_all_IC(:,2), perc)  );
IC_Rampl_std_mean = Rampl_ratio_std_mean_all_IC(idx_score_RR,1);

idx_skew = find(Rampl_skew_all_IC(:,2) < prctile(Rampl_skew_all_IC(:,2), perc)  );
IC_Rampl_skew = Rampl_skew_all_IC(idx_skew,1);

idx_score_kurt = find(Rampl_kurt_all_IC(:,2) < prctile(Rampl_kurt_all_IC(:,2), perc)  );
IC_Rampl_kurt = Rampl_kurt_all_IC(idx_score_kurt,1);

idx_chi2stat = find(Rampl_chi2stat_all_IC(:,2) < prctile(Rampl_chi2stat_all_IC(:,2), perc)  );
IC_Rampl_chi2stat = Rampl_chi2stat_all_IC(idx_chi2stat,1);


% Amplitude of Q peaks
idx_score_Q = find(Qampl_ratio_std_mean_all_IC(:,2) > prctile(Qampl_ratio_std_mean_all_IC(:,2), 100-perc)  );  % ">" and "100-perc" bc heart IC will have a high score
IC_Qampl_std_mean = Qampl_ratio_std_mean_all_IC(idx_score_Q,1);

idx_skew = find(Qampl_skew_all_IC(:,2) > prctile(Qampl_skew_all_IC(:,2), 100-perc)  );
IC_Qampl_skew = Qampl_skew_all_IC(idx_skew,1);

idx_score_kurt = find(Qampl_kurt_all_IC(:,2) < prctile(Qampl_kurt_all_IC(:,2), perc)  );
IC_Qampl_kurt = Qampl_kurt_all_IC(idx_score_kurt,1);

idx_chi2stat = find(Qampl_chi2stat_all_IC(:,2) < prctile(Qampl_chi2stat_all_IC(:,2), perc)  );
IC_Qampl_chi2stat = Qampl_chi2stat_all_IC(idx_chi2stat,1);



% Amplitude of S peaks
idx_score_S = find(Sampl_ratio_std_mean_all_IC(:,2) > prctile(Sampl_ratio_std_mean_all_IC(:,2), 100-perc)  );  % ">" and "100-perc" bc heart IC will have a high score
IC_Sampl_std_mean = Sampl_ratio_std_mean_all_IC(idx_score_S,1);

idx_skew = find(Sampl_skew_all_IC(:,2) > prctile(Sampl_skew_all_IC(:,2), 100-perc)  );
IC_Sampl_skew = Sampl_skew_all_IC(idx_skew,1);

idx_score_kurt = find(Sampl_kurt_all_IC(:,2) < prctile(Sampl_kurt_all_IC(:,2), perc)  );
IC_Sampl_kurt = Sampl_kurt_all_IC(idx_score_kurt,1);

idx_chi2stat = find(Sampl_chi2stat_all_IC(:,2) < prctile(Sampl_chi2stat_all_IC(:,2), perc)  );
IC_Sampl_chi2stat = Sampl_chi2stat_all_IC(idx_chi2stat,1);



% 3) Find IC that meet bpm and SignalAmpl range conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bpm
idx_bpm = intersect(  find((bpm_all_IC(:,2) > bpm_min)),  find(bpm_all_IC(:,2) < bpm_max) );
IC_bpm_ok = bpm_all_IC(idx_bpm,1);

% SignalAmpl range
idx_SignalAmpl_range = find(SignalAmpl_range_all_IC(:,2) < threshold_regularity_signal_minmax);
IC_SignalAmpl_range_ok = SignalAmpl_range_all_IC(idx_SignalAmpl_range,1);



%% Compute final score and find IC with score above threshold (method 1)

%%%%%%%%%%%%%%%%%%%%%%% USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% // To modifiy if you change the nb of conditions below \\
nbr_condition = 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IC_heart_method1 = [];
score_other_IC = [];
all_score = [];
for comp_iter = 1:length(comp.label)

    % Case IC_bug
    if ismember(comp_iter, IC_bug)
        score_other_IC = [score_other_IC, NaN];
        all_score = [all_score; comp_iter, NaN];
    else

        % Check if the IC iter meet the different conditions
            % RR interval conditions
        cond1 = ismember(comp_iter, IC_RR_std_mean);
        cond2 = ismember(comp_iter, IC_RR_skew);
        cond3 = ismember(comp_iter, IC_RR_kurt);
        cond4 = ismember(comp_iter, IC_RR_chi2stat);
            % Rampl conditions
        cond5 = ismember(comp_iter, IC_Rampl_std_mean);
        cond6 = ismember(comp_iter, IC_Rampl_skew);
        cond7 = ismember(comp_iter, IC_Rampl_kurt);
        cond8 = ismember(comp_iter, IC_Rampl_chi2stat);
            % Q and S ampl conditions
        cond9 = ismember(comp_iter, IC_Qampl_std_mean) || ismember(comp_iter, IC_Sampl_std_mean);
        cond10 = ismember(comp_iter, IC_Qampl_skew) || ismember(comp_iter, IC_Sampl_skew);
        cond11 = ismember(comp_iter, IC_Qampl_kurt) || ismember(comp_iter, IC_Sampl_kurt);
        cond12 = ismember(comp_iter, IC_Qampl_chi2stat) || ismember(comp_iter, IC_Sampl_chi2stat);

        % Proportion of conditions met
        prop_cond_ok = (cond1 + cond2 + cond3 + cond4 + cond5 + cond6 + cond7 + cond8 + cond9+cond10+cond11+cond12 ) / nbr_condition ;

        all_score = [all_score; comp_iter, prop_cond_ok];

        % Specific conditions that heart IC MUST meet (bpm and SignalAmplRange)
        cond_bpm = ismember(comp_iter, IC_bpm_ok);
        cond_SignalAmpl_range = ismember(comp_iter, IC_SignalAmpl_range_ok);

        % Determine if it's a heart IC
        if (prop_cond_ok >= threshold_cond_IC_method1) && (cond_bpm == 1) && (cond_SignalAmpl_range == 1)
            IC_heart_method1 = [IC_heart_method1; comp_iter, prop_cond_ok];
            score_other_IC = [score_other_IC, NaN];
        else
            score_other_IC = [score_other_IC, prop_cond_ok];
        end
    end
end

% Extract IC_num, bpm and score for the heart IC
if isempty(IC_heart_method1) == 0
    IC_heart_method1_num = IC_heart_method1(:,1);
    IC_heart_method1_score = IC_heart_method1(:,2);
    IC_heart_method1_bpm = bpm_all_IC(ismember(bpm_all_IC(:,1), IC_heart_method1_num), 2);
else
    IC_heart_method1_num = [];
    IC_heart_method1_score = [];
    IC_heart_method1_bpm = [];
end

% To have an idea of the proportion of conditions met by non-heart IC
Max_score_other_IC = max(score_other_IC);

% For sanity check
Nbr_sec_ICA = size(comp.trial,2);


%% Method 2: find heart IC based on the distribution of scores

% Candidates = IC with score > mean + n*std
IC_heart_method2_candidates = all_score(find( all_score(:,2) > nanmean(all_score(:,2)) + threshold_std_method2*nanstd(all_score(:,2))),1);

% For each candidate IC, check if 1) bpm physiological and 2) the SignalAmplRange
IC_heart_method2 = [];
for i = 1:length(IC_heart_method2_candidates)
    iter = IC_heart_method2_candidates(i);
    if ismember(iter, IC_bpm_ok) && ismember(iter, IC_SignalAmpl_range_ok) && ismember(iter,IC_heart_method2_candidates)
        IC_heart_method2 = [IC_heart_method2; iter];
    end
end


%% Summarize info in a table

% Create a table for recording iter
table_cardiac_IC = table({IC_bug}, {bpm_all_IC(:,2)}, {IC_bpm_ok}, {SignalAmpl_range_all_IC(:,2)}, {IC_SignalAmpl_range_ok},...
    {IC_RR_std_mean}, {IC_RR_skew}, {IC_RR_kurt}, {IC_RR_chi2stat},  ...
    {IC_Rampl_std_mean}, {IC_Rampl_skew}, {IC_Rampl_kurt}, {IC_Rampl_chi2stat}, ...
    {IC_heart_method2}, {IC_heart_method1_num}, {IC_heart_method1_score}, {all_score}, {score_other_IC}, {Max_score_other_IC}, {IC_heart_method1_bpm}, {Nbr_sec_ICA});

% Rename table columns
 col_names = {'IC_bug', 'IC_bpm_all', 'IC_bpm_ok', 'SignalAmpl_range_all', 'SignalAmpl_range_ok',...
        'IC_RR_std_mean', 'IC_RR_skew', 'IC_RR_kurt', 'IC_RR_chi2stat', ...
        'IC_Rampl_std_mean', 'IC_Rampl_skew', 'IC_Rampl_kurt', 'IC_Rampl_chi2stat', ...
        'IC_heart_method2', 'IC_heart_method1', 'Prop_cond_ok_heart_IC', 'Prop_cond_ok_all_IC', 'Prop_cond_ok_other_IC', 'Max_prop_cond_ok_other_IC', 'bpm_final_IC', 'Nbr_sec_ICA'};

table_cardiac_IC.Properties.VariableNames = col_names;



%% Export output


% Extract the cardiac IC depending on the chosing method 
if strcmp(method_chosen, 'absolute_threshold')
    heart_IC = table_cardiac_IC.IC_heart_method1{1}';
elseif strcmp(method_chosen, 'mean_std')
    heart_IC = table_cardiac_IC.IC_heart_method2{1}';
else
    error("Error: method_chosen must take the value 'absolute_threshold' or 'mean_std'")
end



% Plot time course of heart IC

if plot_heart_IC == 1

    % 1) Plot DISTRIBUTION PROP CONDITION OK
    threshold_heart1 = threshold_cond_IC_method1;
    threshold_heart2 = nanmean(all_score(:,2)) + threshold_std_method2*nanstd(all_score(:,2));
    figure
    plot(all_score(:,2),'.','linestyle','none', 'MarkerSize', 15);
    hold on
    xlim([0 length(all_score(:,2))+1]);
    xl = xlim;
    % yl = ylim;
    set(gca,'ylimMode','manual');
    plot(xl, [threshold_heart1, threshold_heart1], 'r--', 'LineWidth', 1);
    plot(xl, [threshold_heart2, threshold_heart2], 'b--', 'LineWidth', 1);
    % Emphasize the rejected IC
    hold on
    emphasizeX = [heart_IC];    % X-coordinates of points to emphasize
    emphasizeY = all_score(heart_IC,2);  % Corresponding Y-coordinates
    plot(emphasizeX, emphasizeY, '.','linestyle','none', 'MarkerSize', 25, 'Color', 'r');
    title(['Prop conditions ok (red: absolute thresh, blue:mean score*' num2str(threshold_std_method2) 'std )'])
    ylabel('Heart score (using pop_iclabel)');
    xlabel('Components');
    % %%%%% SAVE fig as .png %%%%%%%%
    filename = strcat(path_output, '/', file_info, '_ECG_IC_label_fct_find_cardiac.png');
    saveas(gcf, filename);
    close(gcf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 2) If cardiac IC, plot timecourse
    for comp_iter = 1:length(heart_IC)
        figure
        IC_to_plot = heart_IC(comp_iter);
        figure
        plot([1:size(IC_timecourse,2)], IC_timecourse(IC_to_plot,:))
        title(['Cardiac component (IC ' num2str(IC_to_plot) ')'])
        % %%%%% SAVE fig as .png %%%%%%%%
        filename = strcat(path_output, '/', file_info, '_ECG_fct_find_cardiac_time_course_heart_IC ', num2str(IC_to_plot), '.png');
        saveas(gcf, filename);
        close(gcf);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end





end


