function [heart_IC, table_cardiac_IC, aaa_parameters_find_heart_IC] = fct_find_cardiac_IC(comp, method_chosen, plot_heart_IC, path_output, file_info)

%% DOCSTRING

% Author
% % % % % % % % % % % % % % % % % 
% Pierre Champetier (2024)
% Contact: pi.champetier@gmail.com

% Aim
% % % % % % % % % % % % % % % % %
% This function takes EEG independent components (IC) as inputs, and
% identify the cardiac IC. 
% --> To do so!
% 1) The script detects cardiac events (PQRST) in all IC.
% 2) It computes several features (std/mean ratio, skewness, kurtosis...) for the 
% distribution of RR intervals, R, Q and S amplitudes. 
% 3) It identifies the top 3 IC with the lowest values for each feature (or highest
% values depending on the feature), and gives a score (1 or 0) for the feature.
% 4) It computes a global score (sum of score/nb of feature).
% 5) Cardiac IC are identified based on the global score:
%   -method 1: cardiac IC if global score >= threshold (ex: 0.5)
%   -method 2: cardiac if global score >= mean(all_global scores) + n*std(all_global_scores)


% Usage
% % % % % % % % % % % % % % % % % 
% [heart_IC, table_cardiac_IC, parameters_find_heart_IC] = fct_find_cardiac_IC(comp, method_chosen, plot_heart_IC)


% Inputs
% % % % % % % % % % % % % % % % % 
% 1) comp --> output of EEGLAB pop_runica function (i.e., output of the ICA)
% 2) method_chosen --> 'absolute_threshold' (method 1) or 'mean_std' (method 2)
% 3) plot_heart_IC --> 1 or 0 (to plot the IC labelled as cardiac)

% Detail for 'method_chosen'
% % % % % % % % % % % % % % % % % 
% The cardiac IC will be identified 

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



%% Define parameters --> CAN BE CHANGE BY USER

nb_IC_wanted = 3; % The heart IC must be in the top 3 of all IC for skewness, kurtosis... of Rampl, RRintervals...
bpm_min = 45;
bpm_max = 90;
threshold_cond_IC_method1 = 0.5;
threshold_std_method2 = 2.5;
threshold_regularity_signal_minmax = 1.5;
min_recording_duration_sec = 20;

%% Extract parameters

aaa_parameters_find_heart_IC = [];
aaa_parameters_find_heart_IC.nb_IC_wanted = nb_IC_wanted;
aaa_parameters_find_heart_IC.bpm_min = bpm_min;
aaa_parameters_find_heart_IC.bpm_max = bpm_max;
aaa_parameters_find_heart_IC.threshold_regularity_signal_minmax = threshold_regularity_signal_minmax;
aaa_parameters_find_heart_IC.threshold_cond_IC_method1 = threshold_cond_IC_method1;
aaa_parameters_find_heart_IC.threshold_std_method2 = threshold_std_method2;
aaa_parameters_find_heart_IC.min_recording_duration_sec = min_recording_duration_sec;

fs = comp.fsample;


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

    nbr_mini_segments = 10;

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
IC_bpm = bpm_all_IC(idx_bpm,1);

% SignalAmpl range
idx_SignalAmpl_range = find(SignalAmpl_range_all_IC(:,2) < threshold_regularity_signal_minmax);
IC_SignalAmpl_range = SignalAmpl_range_all_IC(idx_SignalAmpl_range,1);



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
        cond_bpm = ismember(comp_iter, IC_bpm);
        cond_SignalAmpl_range = ismember(comp_iter, IC_SignalAmpl_range);

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

% For each candidate IC, check if 1) bpm physiological and 2) the
IC_heart_method2 = [];
for i = 1:length(IC_heart_method2_candidates)
    iter = IC_heart_method2_candidates(i);
    if ismember(iter, IC_bpm) && ismember(iter, IC_SignalAmpl_range) && ismember(iter,IC_heart_method2_candidates)
        IC_heart_method2 = [IC_heart_method2; iter];
    end
end


%% Summarize info in a table

% Create a table for recording iter
table_cardiac_IC = table({IC_bug}, {IC_bpm}, ...
    {IC_RR_std_mean}, {IC_RR_skew}, {IC_RR_kurt}, {IC_RR_chi2stat},  ...
    {IC_Rampl_std_mean}, {IC_Rampl_skew}, {IC_Rampl_kurt}, {IC_Rampl_chi2stat}, ...
    {IC_heart_method2}, {IC_heart_method1_num}, {IC_heart_method1_score}, {all_score}, {score_other_IC}, {Max_score_other_IC}, {IC_heart_method1_bpm}, {Nbr_sec_ICA});

% Rename table columns
 col_names = {'IC_bug', 'IC_bpm', ...
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
    filename = strcat(path_output, '/', file_info, '_ECG_IC_label_PIERRE_fct.png');
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
        filename = strcat(path_output, '/', file_info, '_ECG_PIERRE_time_course_heart_IC ', num2str(IC_to_plot), '.png');
        saveas(gcf, filename);
        close(gcf);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end





end


