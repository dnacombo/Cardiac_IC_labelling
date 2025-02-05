function [aaa_parameters_find_heart_IC, output_for_zscore_corMatrix_ROC, output_for_user] = CARACAS(cfg, comp)

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
%   -method_chosen --> 'absolute_threshold' (method 1) or 'mean_std' (method 2)
%   -plot_heart_IC --> 1 or 0 (to plot the IC labelled as cardiac)
%   -path_output --> path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
%   -file_info --> name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
%   -nb_IC_wanted --> number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
%   -bpm_min and bpm_max  --> expected heart beat per min, for sanity check [default: 45 and 90]
%   -threshold_cond_IC_method1 --> minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
%   -threshold_std_method2 --> if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
%   -min_recording_duration_sec --> minimum duration (in sec) of the IC timecourse (default: 20]
%   -mini_bouts_duration_for_SignalAmplRange --> for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
% - threshold_regularity_signal_minmax --> For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]
% 
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
% A_fct_test_unif_NEW.m (To test the uniform distribution of cardiac events)
% ECG_PQRST_VERSION_3_PC (Toolbox used to detect cardiac events) --> Rq: I did small modifications of 'compute_fudicial_peaks_live17_c' and 'preprocess_window_ecg' functions to avoid errors with some IC (all modifications are preceeded by a comment "PIERRE").

% add ECG toolbox to the path
% addpath(fullfile(fileparts(which(mfilename)),'ECG_PQRST_VERSION_3_PC'))
addpath(fullfile(fileparts(which(mfilename)),'heart_functions'))
%% CONSTANT PARAMETERS

fs = comp.fsample;

% Window for ERP
window = 0.2*fs; % 0.2s before and 0.2s after the R peak

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
    bpm_min = 35;
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
    threshold_cond_IC_method1 = 0.6;
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

% IC to not analyze (because user has already labelled it for instance)
if isfield(cfg, 'IC_to_not_analyze') == 1
    IC_to_not_analyze = cfg.IC_to_not_analyze;
    % Check user has given number(s) for IC to not analyze
    if isnumeric(IC_to_not_analyze) == 0 && ~isequal(IC_to_not_analyze, [])
        error('Error in cfg.IC_to_not_analyze: please enter a number (or use [] /do not use cfg.IC_to_not_analyze).')
    end
    % Check user has given an IC that exists
    if sum(ismember(IC_to_not_analyze, [1:length(comp.label)])) ~= length(IC_to_not_analyze)
        error('Error in cfg.IC_to_not_analyze: you want to remove IC that does not exist.')
    end
else
    IC_to_not_analyze = []; 
end

%% Pre-defined col_names (as it might be needed in "Extract parameters section if recording_duration_sec < min_recording_duration_sec)

col_names = {'IC_bug', 'IC_bpm_all', 'IC_bpm_ok', 'SignalAmpl_range_all', 'SignalAmpl_range_ok',...
        'IC_ERPampl_median',  'IC_ERPampl_std', 'IC_ERPampl_skew', 'IC_ERPampl_kurt', 'IC_ERPampl_chi2stat', 'IC_ERPampl_std_median',...
        'IC_RR_std_mean', 'IC_RR_skew', 'IC_RR_kurt', 'IC_RR_chi2stat', ...
        'IC_Rampl_std_mean', 'IC_Rampl_skew', 'IC_Rampl_kurt', 'IC_Rampl_chi2stat', ...
        'IC_heart_method2', 'IC_heart_method1', 'Prop_cond_ok_heart_IC', 'Prop_cond_ok_all_IC', 'Prop_cond_ok_other_IC', 'Max_prop_cond_ok_other_IC', 'bpm_final_IC', 'Nbr_sec_ICA'};


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
aaa_parameters_find_heart_IC.IC_to_not_analyze = IC_to_not_analyze;

% Checking that recording_duration_sec > mini_bouts_duration_for_SignalAmplRange
sample_tot = 0;
for i = 1:length(comp.trial)
    sample_tot = sample_tot + length(comp.trial{1,i});
end
recording_duration_sec = sample_tot / fs;

if recording_duration_sec < mini_bouts_duration_for_SignalAmplRange
    error('Error: attempt to divide the recording into bouts with longer than the original recording (for the Signal Ampl Range checking). Change the value of mini_bouts_duration_for_SignalAmplRange variable in the function.')
end

% Checking that recording_duration_sec > min_recording_duration_sec
if recording_duration_sec < min_recording_duration_sec
    warning('recording_duration_sec < min_recording_duration_sec')
    table_cardiac_IC = table({[1:length(comp.label)]}, {[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]});
    table_cardiac_IC.Properties.VariableNames = col_names;
    heart_IC = [];
else



    %% OLD: Extract timecourse of each IC for all 1s-mini trials

    IC_timecourse = [];
    for i = 1:size(comp.trial,2)
        IC_timecourse = [IC_timecourse, comp.trial{1,i}];
    end

    %% NEW: Recreate continuous comp data
    cfgtmp = [];
    if isfield(comp, 'sampleinfo') 
        cfgtmp.trl = [1, comp.sampleinfo(end,2), 0];
    else
        % error('Error: no sampleinfo in comp')
        cfgtmp.trl = [1, sum(cellfun(@numel,comp.time)), 0]; % Start, End, Offset
    end
    comp_continu = ft_redefinetrial(cfgtmp, comp);

    % Extract timecourse
    timecourse_all_IC = comp_continu.trial{1};

    %% Plot IC
    % for comp_iter = 1:6
    %     figure
    %     plot([1:length(IC_timecourse)], IC_timecourse(comp_iter,:))
    %     title(num2str(comp_iter))
    % end
    % % 


    IC_timecourse_new = comp_continu.trial{1};
    %
    % for comp_iter = 1:6
    %     figure
    %     plot([1:size(IC_timecourse_new,2)], IC_timecourse_new(comp_iter,:))
    %     title(num2str(comp_iter))
    % end


    %% DETECT CARDIAC EVENTS

    plot_data = 0; % (see FUNCTION DOCUMENTATION)

    IC_bug = [];
    IC_with_less_than_10_R_peak = [];

    bpm_all_IC = [];

    % RR features
    RR_mean_all = [];
    RR_std_all = [];
    RR_skew_all = [];
    RR_kurt_all = [];
    RR_ratio_std_mean_all_IC = [];
    RR_chi2stat_all_IC = [];

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

    ERP_average_for_each_IC = [];
    ERP_details_for_each_IC = [];
    ERPampl_median = [];
    ERPampl_std = [];
    ERPampl_std_median = [];
    ERPampl_skew = [];
    ERPampl_kurt = [];
    ERPampl_chi2stat = [];



    nbr_cardiac_event = [];
    IC_not_cardiac_bc_PQstd = [];
    IC_not_cardiac_bc_RRstd = [];
    IC_not_cardiac_bc_Ramplstd = [];

    for comp_iter = 1:length(comp.label)
        % Skip if user don't want to analyze this IC
        if ismember(comp_iter, IC_to_not_analyze)
            ERP_average_for_each_IC = [ERP_average_for_each_IC; repmat(NaN,1,window*2+1);];
            ERP_details_for_each_IC{1,comp_iter} = {};
            ERPampl_median = [ERPampl_median; comp_iter, NaN];

            bpm_all_IC = [bpm_all_IC; comp_iter, NaN];
            SignalAmpl_range_all_IC = [SignalAmpl_range_all_IC; comp_iter, NaN];

            RR_mean_all = [RR_mean_all; comp_iter, NaN];
            RR_std_all = [RR_std_all; comp_iter, NaN];
            
            RR_ratio_std_mean_all_IC = [RR_ratio_std_mean_all_IC; comp_iter, NaN];
            RR_skew_all = [RR_skew_all; comp_iter, NaN];
            RR_kurt_all = [RR_kurt_all; comp_iter, NaN];
            RR_chi2stat_all_IC = [RR_chi2stat_all_IC; comp_iter, NaN];
            Rampl_ratio_std_mean_all_IC = [Rampl_ratio_std_mean_all_IC; comp_iter, NaN];
            Rampl_skew_all_IC = [Rampl_skew_all_IC; comp_iter, NaN];
            Rampl_kurt_all_IC = [Rampl_kurt_all_IC; comp_iter, NaN];
            Rampl_chi2stat_all_IC = [Rampl_chi2stat_all_IC; comp_iter, NaN];
            Qampl_ratio_std_mean_all_IC = [Qampl_ratio_std_mean_all_IC; comp_iter, NaN];
            Qampl_skew_all_IC = [Qampl_skew_all_IC; comp_iter, NaN];
            Qampl_kurt_all_IC = [Qampl_kurt_all_IC; comp_iter, NaN];
            Qampl_chi2stat_all_IC = [Qampl_chi2stat_all_IC; comp_iter, NaN];
            Sampl_ratio_std_mean_all_IC = [Sampl_ratio_std_mean_all_IC; comp_iter, NaN];
            Sampl_skew_all_IC = [Sampl_skew_all_IC; comp_iter, NaN];
            Sampl_kurt_all_IC = [Sampl_kurt_all_IC; comp_iter, NaN];
            Sampl_chi2stat_all_IC = [Sampl_chi2stat_all_IC; comp_iter, NaN];
            continue
        end

        % Select IC iter
        % ecg_1 = IC_timecourse(comp_iter,:);

        % Detect cardiac events
        % try
            % ecg_f = preprocess_window_ecg(ecg_1, fs);
            % ecg_f = ecg_1;
            % [locs_P,locs_Q,locs_R,locs_S,locs_T] = compute_fudicial_peaks_live17_c(ecg_f, fs, plot_data); % Why the nb of samples in ecg_f is lower than in ecg_1?? For some IC it's super low....
           
            cfg_peak = [];
            % cfg_peak.plotall         = 1;
            % cfg_peak.plotthresh      = 0;
            % cfg_peak.plotbeat        = 1;
            % cfg_peak.plotcorr        = 0;
            % cfg_peak.plotfinal       = 1;
            cfg_peak.channel = comp.label{comp_iter};
 
            try

                [HeartBeats] = heart_peak_detect(cfg_peak,comp_continu);


                locs_P = [HeartBeats.P_sample];
                locs_Q = [HeartBeats.Q_sample];
                locs_R = [HeartBeats.R_sample];
                locs_S = [HeartBeats.S_sample];
                locs_T = [HeartBeats.T_sample];


        


                % PQ intervals
                %%%%%%%%%%%%%
                PQ_intervals_sec = (locs_P - locs_Q)/comp_continu.fsample;

                % figure;
                % hist(PQ_intervals_sec)

                low_threshold = prctile(PQ_intervals_sec, 15);  % 5th percentile
                high_threshold = prctile(PQ_intervals_sec, 85); % 95th percentile

                filtered_PQ_intervals = PQ_intervals_sec(PQ_intervals_sec >= low_threshold & PQ_intervals_sec <= high_threshold);

                % figure;
                % hist(filtered_PQ_intervals);

                if std(filtered_PQ_intervals) > abs(mean(filtered_PQ_intervals))/3
                    IC_not_cardiac_bc_PQstd = [IC_not_cardiac_bc_PQstd, comp_iter];
                end



                % RR intervals
                %%%%%%%%%%%%%%
                RR_intervals_sec = diff(locs_R)/comp_continu.fsample;
                % figure;
                % hist(RR_intervals_sec);

                % Remove lowest 5% and highest 20% lowest values (that biase the std)
                low_threshold = prctile(RR_intervals_sec, 0);  
                high_threshold = prctile(RR_intervals_sec, 70); 

                filtered_RR_intervals = RR_intervals_sec(RR_intervals_sec >= low_threshold & RR_intervals_sec <= high_threshold);
                
                % figure;
                % hist(filtered_RR_intervals);

                if std(filtered_RR_intervals) > mean(filtered_RR_intervals)/3
                    IC_not_cardiac_bc_RRstd = [IC_not_cardiac_bc_RRstd, comp_iter];
                end


                % Rampl
                %%%%%%%
                time_course_IC_iter = timecourse_all_IC(comp_iter,:);
                Rampl = time_course_IC_iter(locs_R);
                Rampl = abs(Rampl);

                % figure;
                % hist(Rampl);

                % Remove lowest 15% and highest 15% lowest values (that biase the std)
                low_threshold = prctile(Rampl, 15);  % 5th percentile
                high_threshold = prctile(Rampl, 85); % 95th percentile

                filtered_Rampl = Rampl(Rampl >= low_threshold & Rampl <= high_threshold);

                % figure;
                % hist(filtered_Rampl);

                if std(filtered_Rampl) > abs(mean(filtered_Rampl))/3
                    IC_not_cardiac_bc_Ramplstd = [IC_not_cardiac_bc_Ramplstd, comp_iter];
                end                



            catch
                warning('error in heart_peak_detect')
                % nbr_cardiac_event = [nbr_cardiac_event; comp_iter 0];
                locs_P = [];
            end

            nbr_cardiac_event = [nbr_cardiac_event; comp_iter length(locs_P)];




            %%  Metric to check homogeneous SignalAmpl across the recording
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            % Remove NaN
            time_course_comp_iter = timecourse_all_IC(comp_iter,:);
            time_course_comp_iter(find(isnan(time_course_comp_iter))) = [];


            % Divide the recordings into mini segments and check that the abs(mean) is similar for each segment
            nbr_mini_segments = round((length(time_course_comp_iter) / fs) / mini_bouts_duration_for_SignalAmplRange);

            % Extract the absolute value of the signal
            array = abs(time_course_comp_iter);

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


    % 1) bpm
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % Option 1: with data-driven threshold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold_IC_4std = median(nbr_cardiac_event(:,2)) +  4*std(nbr_cardiac_event(:,2));
    threshold_IC_3std = median(nbr_cardiac_event(:,2)) +  3*std(nbr_cardiac_event(:,2));


    threshold_IC = threshold_IC_3std;
    % heart_IC = nbr_cardiac_event(find(nbr_cardiac_event(:,2) > threshold_IC),1);


    % Option 2: with physiological bpm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbr_cardiac_event(:,3) = nbr_cardiac_event(:,2)/(recording_duration_sec /60);
    heart_IC = nbr_cardiac_event(nbr_cardiac_event(:,3) >= bpm_min & nbr_cardiac_event(:,3) <= bpm_max, :);
    heart_IC = heart_IC(:,1);



    % 2) Remove potential cardiac if too high std for RR interval or Rampl
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    % For PQ interval
    idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_PQstd));
    heart_IC(idx_IC_with_too_high_std) = [];    

    % % For QR interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_QRstd));
    % heart_IC(idx_IC_with_too_high_std) = [];
    % 
    % % For RS interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_RSstd));
    % heart_IC(idx_IC_with_too_high_std) = [];
    % 
    % % For ST interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_STstd));
    % heart_IC(idx_IC_with_too_high_std) = [];

    % For RR interval
    idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_RRstd));
    heart_IC(idx_IC_with_too_high_std) = [];

        % For Rampl
    idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_Ramplstd));
    heart_IC(idx_IC_with_too_high_std) = [];


    % 3) Remove potential cardiac if too high SignalAmpl range
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    idx_SignalAmpl_range = find(SignalAmpl_range_all_IC(:,2) > threshold_regularity_signal_minmax);
    IC_not_cardiac_bc_Ramplstd = SignalAmpl_range_all_IC(idx_SignalAmpl_range,1);

    idx_IC_with_too_high_SignalAmpl_range = find(ismember(heart_IC, IC_not_cardiac_bc_Ramplstd));
    heart_IC(idx_IC_with_too_high_SignalAmpl_range) = [];


    %% Save output
    
    output_for_zscore_corMatrix_ROC =[];

    output_for_user = [];
    output_for_user.threshold_IC_3std = threshold_IC_3std;
    output_for_user.threshold_IC_4std = threshold_IC_4std;

    output_for_user.threshold_IC = threshold_IC;
    output_for_user.IC_nbr_cardiac_event = nbr_cardiac_event(:,2);

    if isempty(heart_IC) == 0
        output_for_user.heart_IC = heart_IC;
        output_for_user.heart_IC_nbr_cardiac_event = nbr_cardiac_event(find(nbr_cardiac_event(:,2) > threshold_IC),2);
    else
        output_for_user.heart_IC = NaN;
        output_for_user.heart_IC_nbr_cardiac_event = NaN;
    end

    %%
        % catch
        %     IC_bug = [IC_bug, comp_iter]; % I don't know why but the function preprocess_window_ecg bug for some recordings
        %     continue;
        % end


%         % Check enough data
%         duration_sec_data_iter = length(ecg_1) / fs;
%         if duration_sec_data_iter < min_recording_duration_sec
%             IC_bug = [IC_bug, comp_iter];
%             continue
%         end
% 
% 
%         %% Plot cardiac events
% 
%         % figure(1);
%         % hold on
%         % plot(ecg_f,'linewidth',1.5);
%         % plot(locs_R,ecg_f(locs_R),'rv','MarkerFaceColor','r');
%         % plot(locs_Q,ecg_f(locs_Q),'rs','MarkerFaceColor','g');
%         % plot(locs_S,ecg_f(locs_S),'o','MarkerFaceColor','y');
%         % plot(locs_P,ecg_f(locs_P),'v','MarkerFaceColor','m');
%         % plot(locs_T,ecg_f(locs_T),'s','MarkerFaceColor','k');
%         % % plot(T2,'linewidth',2,'Linestyle','--','Color','red');
%         % % plot(T1,'linewidth',2,'Linestyle','--','Color','green');
%         % legend('ECG Signal','R peaks','Q Peaks','S Peaks','P peaks','T peaks','T-P1','T-P2');
%         % hold off
%         % title('QRS complex being marked on the fitered ECG signal along with the P and T peaks');
%         % xlabel('Samples');
%         % ylabel('Normalised Amplitude');
%         % grid on
% 
% 
%         %% ERP cardiac events (QRS peak-to-peak amplitude)
% 
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%         % 1) Choose data :
%         %%%%%%%%%%%%%%%%%%
%             % Option 1: ecg_1 (original IC_timecourse) 
%             % Option 2: ecg_f (filtered output of preprocess_window_ecg)
%         data_for_ERP = ecg_f; 
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%         % 
%         % figure; plot(data_for_ERP(locs_R(R_iter)-window : locs_R(R_iter)+window))
%         % figure; plot(ecg_1(locs_R(R_iter)-10*window : locs_R(R_iter)+10*window));
%         % figure; plot(ecg_f(locs_R(R_iter)-10*window : locs_R(R_iter)+10*window))
% 
% 
%         % 2) Extract signal around all R peaks
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         all_cardiac_events = [];
%         cardiac_events_iter = [];
% 
%         for R_iter = 1:length(locs_R)
%             % Try in case the first R is at the very begining of at the very end of the recording (so we cannot extract the signal before or after)
%             try
%                 cardiac_events_iter = [cardiac_events_iter; data_for_ERP(locs_R(R_iter)-window : locs_R(R_iter)+window)];
%             catch
%                 continue
%             end
%         end
% 
%         ERP_details_for_each_IC{1,comp_iter} = cardiac_events_iter;
% 
% 
%         if size(cardiac_events_iter,1) > 10
%             % Compute the ERP (median of cardiac events)
%             ERP_median_ampl_iter = median(cardiac_events_iter); 
% 
%             % 3) Compute peak-to-peak ERP_ampl and gather with other IC
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % figure; plot(ERP_median_ampl_iter)
% 
%             % 3.1) Find sample of positive peak (R peak)
%             idx_pos_peak = round(size(ERP_median_ampl_iter,2)/2); % = middle of the window
% 
%             % 3.2) Find sample of negative peak
%             % 3.2.1) Define begin and end of the searching window to find the negative peal lowest value (negative peak)
%             beginning = round(idx_pos_peak-window/3);
%             ending = round(idx_pos_peak+window/3);
%             % 3.2.2) Find the negative peak within the window
%             [minValue, minIdx] = min(ERP_median_ampl_iter(beginning:ending));
%             idx_neg_peak = beginning+minIdx-1;
% 
%             % 3.3) Compute the peak-to-peak amplitude of each QRS complex
%             ERP_P2P_iter = cardiac_events_iter(:,idx_pos_peak) - cardiac_events_iter(:,idx_neg_peak);
% 
%             % 3.4) Compute metrics (median, std, skewness, kurtosis) of the QRS amplitude
%             ERPampl_median_iter = median(ERP_P2P_iter);
%             ERPampl_std_iter = std(ERP_P2P_iter);
%             ERPampl_skew_iter = skewness(ERP_P2P_iter);
%             ERPampl_kurt_iter = kurtosis(ERP_P2P_iter);
% 
% 
%         % Case not enough cardiac event
%         else
%             ERP_median_ampl_iter = repmat(NaN,1,window*2+1);
%             % ERP_std_ampl_iter = repmat(NaN,1,window*2+1);
%             ERPampl_median_iter = NaN;
%             ERPampl_std_iter = NaN;
%             ERPampl_skew_iter = NaN;
%             ERPampl_kurt_iter = NaN;
%         end
% 
% 
%         % 4) Gather with other subjects
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ERP_average_for_each_IC = [ERP_average_for_each_IC; ERP_median_ampl_iter];
% 
%         ERPampl_median = [ERPampl_median; comp_iter, ERPampl_median_iter];
%         ERPampl_std = [ERPampl_std; comp_iter, ERPampl_std_iter];
%         ERPampl_std_median = [ERPampl_std_median; comp_iter, ERPampl_std_iter/ERPampl_median_iter];
%         ERPampl_skew = [ERPampl_skew; comp_iter, ERPampl_skew_iter];
%         ERPampl_kurt = [ERPampl_kurt; comp_iter, ERPampl_kurt_iter];
% 
%         % 5) Test uniform distribution of RR intervals
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if size(cardiac_events_iter,1) > 10  
%             chi2stat = A_fct_test_unif_NEW(ERP_P2P_iter);
%         else
%             chi2stat = NaN;
%         end
%         ERPampl_chi2stat = [ERPampl_chi2stat; comp_iter, chi2stat];
% 
% 
%         % 6) Plot std and median ERP
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         % figure;
%         % IC_to_plot = 3;
%         % for i = 1:size(ERP_details_for_each_IC{1, IC_to_plot}, 1)
%         %     hold on
%         %     plot(ERP_details_for_each_IC{1,IC_to_plot}(i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
%         % end
%         % plot(ERP_median_ampl_iter, 'Color', 'green', 'LineWidth', 2)
%         % plot(ERP_std_ampl_iter, 'Color', 'blue', 'LineWidth', 2)
%         % plot(ERP_std_ampl_iter./ERP_median_ampl_iter, 'Color', 'red', 'LineWidth', 2)
% 
% 
%         %% Characteristics of RR intervals
% 
%         % 1) std/mean ratio, skew, kurt
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % Computation
%         % % % % % % %
% 
%         % All parameters available --> could be also used to improve the nb of conditions to check and then improve performance?
%         % [mean_QR, mean_SS, mean_ST, mean_RS, std_QR, std_SS, std_ST, skew_QR, skew_SS, skew_ST, skew_RS, kurt_QR, kurt_SS, kurt_ST, kurt_RS, std_RS] = compute_mean_intervals2_c(locs_R,locs_S,locs_Q,locs_T,fs);
%         % [RR_mean, QRS_mean, mean_QT, mean_PR, RR_std, RR_skew, RR_kurt, QRS_std, QRS_skew, QRS_kurt, std_PR, skew_PR, kurt_PR, std_QT, skew_QT, kurt_QT]= compute_mean_interval_c(locs_P,locs_Q,locs_R,locs_S,locs_T,fs);
% 
%         % Those we need
%         % OLD WRONG: [RR_mean, RR_std, RR_skew, RR_kurt]= compute_mean_interval_c(locs_P,locs_Q,locs_R,locs_S,locs_T,fs);
%         [RR_mean,QRS_mean,mean_QT,mean_PR,RR_std,RR_skew,RR_kurt,QRS_std,QRS_skew,QRS_kurt,std_PR,skew_PR,kurt_PR,std_QT,skew_QT,kurt_QT]= compute_mean_interval_c(locs_P,locs_Q,locs_R,locs_S,locs_T,fs);
% 
%         nbr_cardiac_event = length(locs_R);
%         bpm = 60 * nbr_cardiac_event/(duration_sec_data_iter);
% 
%         % Extract RR interval
%         RR_interval = diff(locs_R) /fs;
% 
%         % Gather with other IC
%         % % % % % % % % % % % %
%         RR_mean_all = [RR_mean_all; comp_iter, RR_mean];
%         RR_std_all = [RR_std_all; comp_iter, RR_std];
%         RR_ratio_std_mean_all_IC = [RR_ratio_std_mean_all_IC; comp_iter, RR_std/RR_mean];
%                 % RR_ratio_std_mean_all_IC = [RR_ratio_std_mean_all_IC; comp_iter, RR_std/median(RR_interval)] %
% 
%         RR_skew_all = [RR_skew_all; comp_iter, RR_skew];
%         RR_kurt_all = [RR_kurt_all; comp_iter, RR_kurt];
%         bpm_all_IC = [bpm_all_IC; comp_iter, bpm];
% 
% 
% 
%         % 2) Test uniform distribution of RR intervals
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%         if length(locs_R) > 10  
%             chi2stat = A_fct_test_unif_NEW(locs_R);
%         else
%             chi2stat = NaN;
%         end
%         RR_chi2stat_all_IC = [RR_chi2stat_all_IC; comp_iter, chi2stat];
% 
%         %% Plot for checking distribution of ERPampl and RR_interval
% 
%         % % ERP_ampl
%         % figure;
%         % histogram(ERP_P2P_iter, 'Normalization', 'probability');
%         % title(['Comp ' num2str(comp_iter) ' (std = '  num2str(std(ERP_P2P_iter))  ', skew='  num2str(skewness(ERP_P2P_iter)) ')']); 
%         % xlabel('ERP ampl');
%         % ylabel('Proportion')
% 
%         % % RR_interval
%         % figure; 
%         % histogram(RR_interval)%, 'Normalization', 'probability');
%         % title(['Comp ' num2str(comp_iter) ' (std='  num2str(std(RR_interval)) ', skew='  num2str(skewness(RR_interval)) ')']); 
%         % xlabel('RR interval (in sec)'); 
%         % ylabel('Proportion')
% 
% 
% 
%         %% Characteristics of the AMPLITUDE of R, S and Q peaks
% 
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%         %                            R PEAKS
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%         % 1) std/mean ratio, skew and kurt
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Rampl_ratio_std_mean_all_IC = [Rampl_ratio_std_mean_all_IC; comp_iter,  std(ecg_f(locs_R)) / mean(ecg_f(locs_R))]; 
%         Rampl_kurt_all_IC = [Rampl_kurt_all_IC; comp_iter, kurtosis(ecg_f(locs_R))];
%         Rampl_skew_all_IC = [Rampl_skew_all_IC; comp_iter, skewness(ecg_f(locs_R))];
% 
%         % 2) Test uniform distribution of amplitude of R peaks
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if length(locs_R) > 10  
%             chi2stat = A_fct_test_unif_NEW(ecg_f(locs_R));
%         else
%             chi2stat = NaN;
%         end
%         Rampl_chi2stat_all_IC = [Rampl_chi2stat_all_IC; comp_iter, chi2stat];
% 
% 
% 
% 
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%         %                            Q PEAKS
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%         % 1) std/mean ratio, skew and kurt
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Qampl_ratio_std_mean_all_IC = [Qampl_ratio_std_mean_all_IC; comp_iter, std(ecg_f(locs_Q)) / mean(ecg_f(locs_Q))];
%         Qampl_kurt_all_IC = [Qampl_kurt_all_IC; comp_iter, kurtosis(ecg_f(locs_Q))];
%         Qampl_skew_all_IC = [Qampl_skew_all_IC; comp_iter, skewness(ecg_f(locs_Q))];
% 
%         % 2) Test uniform distribution of amplitude of R peaks
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if length(locs_R) > 10  
%             chi2stat = A_fct_test_unif_NEW(ecg_f(locs_Q));
%         else
%             chi2stat = NaN;
%         end
%         Qampl_chi2stat_all_IC = [Qampl_chi2stat_all_IC; comp_iter, chi2stat];
% 
% 
% 
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%         %                            S PEAKS
%         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
%         % 1) std/mean ratio, skew and kurt
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Sampl_ratio_std_mean_all_IC = [Sampl_ratio_std_mean_all_IC; comp_iter, std(ecg_f(locs_S)) / mean(ecg_f(locs_S))];
%         Sampl_kurt_all_IC = [Sampl_kurt_all_IC; comp_iter, kurtosis(ecg_f(locs_S))];
%         Sampl_skew_all_IC = [Sampl_skew_all_IC; comp_iter, skewness(ecg_f(locs_S))];
% 
%         % 2) Test uniform distribution of amplitude of R peaks
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if length(locs_R) > 10  
%             chi2stat = A_fct_test_unif_NEW(ecg_f(locs_S));
%         else
%             chi2stat = NaN;
%         end
%         Sampl_chi2stat_all_IC = [Sampl_chi2stat_all_IC; comp_iter, chi2stat];
% 
% 
% 
        % %% Metric to check homogeneous SignalAmpl across the recording
        % % Divide the recordings into mini segments and check that the abs(mean)
        % % is similar for each segment
        % 
        % nbr_mini_segments = round((length(IC_timecourse(comp_iter,:)) / fs) / mini_bouts_duration_for_SignalAmplRange);
        % 
        % % Extract the absolute value of the signal
        % array = abs(IC_timecourse(comp_iter,:));
        % 
        % % Determine the size of each small array
        % sizeOfSmallArray = floor(length(array)/ nbr_mini_segments);
        % 
        % % Preallocate cell array for the small arrays
        % smallArrays = cell(1, nbr_mini_segments);
        % 
        % % Fill the small arrays
        % for i = 1:nbr_mini_segments
        %     startIdx = (i-1) * sizeOfSmallArray + 1;
        %     if i < nbr_mini_segments
        %         endIdx = i * sizeOfSmallArray;
        %     else
        %         % Last array takes any remaining elements
        %         endIdx = length(array);
        %     end
        %     smallArrays{i} = array(startIdx:endIdx);
        % end
        % 
        % % Compute mean for each segment
        % SignalAmpl_mean_segments = [];
        % for i = 1:nbr_mini_segments
        %     SignalAmpl_mean_segments = [SignalAmpl_mean_segments, mean(smallArrays{i})];
        % end
        % 
        % % Compute a metric to evaluate if the signal is regular over the recording (max-min)/mean
        % SignalAmpl_range = (max(SignalAmpl_mean_segments) - min(SignalAmpl_mean_segments)) / min(SignalAmpl_mean_segments);
        % SignalAmpl_range_all_IC = [SignalAmpl_range_all_IC; comp_iter, SignalAmpl_range];

% 
%         %% Check that the features were computed on at least 10 cardiac events
% 
%         if nbr_cardiac_event < 10
%             IC_with_less_than_10_R_peak = [IC_with_less_than_10_R_peak, comp_iter];
% 
%             RR_ratio_std_mean_all_IC(find(RR_ratio_std_mean_all_IC(:,1) == comp_iter), 2) = NaN;
%             RR_skew_all(find(RR_skew_all(:,1) == comp_iter), 2) = NaN;            
%             RR_kurt_all(find(RR_kurt_all(:,1) == comp_iter), 2) = NaN;
%             RR_chi2stat_all_IC(find(RR_chi2stat_all_IC(:,1) == comp_iter), 2) = NaN;
%             Rampl_ratio_std_mean_all_IC(find(Rampl_ratio_std_mean_all_IC(:,1) == comp_iter), 2) = NaN;
%             Rampl_skew_all_IC(find(Rampl_skew_all_IC(:,1) == comp_iter), 2) = NaN;
%             Rampl_kurt_all_IC(find(Rampl_kurt_all_IC(:,1) == comp_iter), 2) = NaN;
%             Rampl_chi2stat_all_IC(find(Rampl_chi2stat_all_IC(:,1) == comp_iter), 2) = NaN;
%             Qampl_ratio_std_mean_all_IC(find(Qampl_ratio_std_mean_all_IC(:,1) == comp_iter), 2) = NaN;
%             Qampl_skew_all_IC(find(Qampl_skew_all_IC(:,1) == comp_iter), 2) = NaN;
%             Qampl_kurt_all_IC(find(Qampl_kurt_all_IC(:,1) == comp_iter), 2) = NaN;
%             Qampl_chi2stat_all_IC(find(Qampl_chi2stat_all_IC(:,1) == comp_iter), 2) = NaN;
%             Sampl_ratio_std_mean_all_IC(find(Sampl_ratio_std_mean_all_IC(:,1) == comp_iter), 2) = NaN;
%             Sampl_skew_all_IC(find(Sampl_skew_all_IC(:,1) == comp_iter), 2) = NaN;
%             Sampl_kurt_all_IC(find(Sampl_kurt_all_IC(:,1) == comp_iter), 2) = NaN;
%             Sampl_chi2stat_all_IC(find(Sampl_chi2stat_all_IC(:,1) == comp_iter), 2) = NaN;
%         end
% 
% 
%     end
% 
% 
%         %% Homogenize metrics (to have always "heart IC should have low value" for each metric)
% 
%         % Ex: heart IC will have a high raw Qampl_ratio_std_mean_all_IC (bc Q peak is negative, so mean(Qampl)<0, so high std/mean) ==> so we multiply by -1
%         ERPampl_median(:,2) = -ERPampl_median(:,2);
%         ERPampl_std(:,2) = -ERPampl_std(:,2); % High values of std for cardiac IC bc higher values of amplitude
% 
%         Qampl_ratio_std_mean_all_IC(:,2) = -Qampl_ratio_std_mean_all_IC(:,2);
%         Qampl_skew_all_IC(:,2) = -Qampl_skew_all_IC(:,2);
%         Sampl_ratio_std_mean_all_IC(:,2) = -Sampl_ratio_std_mean_all_IC(:,2);
%         Sampl_skew_all_IC(:,2) = -Sampl_skew_all_IC(:,2);
% 
% 
% 
%     %% Compute metrics to find heart IC 
% 
%     % 0) Manual visualization
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     % figure
%     % plot(Rampl_kurt_all_IC(:,1), Rampl_kurt_all_IC(:,2))
%     % figure
%     % plot(RR_ratio_std_mean_all_IC(:,1), RR_ratio_std_mean_all_IC(:,2))
%     % figure
%     % plot(Rampl_skew_all_IC(:,1), Rampl_skew_all_IC(:,2))
% 
%     % [sortedAbsArray, sortedIndices] = sort(abs(skew_RR_all_IC(:,2)));
% 
% 
%     % 1) Add zscores in column 3
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ERPampl_median = [ERPampl_median, nanzscore(ERPampl_median(:,2))];
%     ERPampl_std = [ERPampl_std, nanzscore(ERPampl_std(:,2))];
%     ERPampl_std_median = [ERPampl_std_median, nanzscore(ERPampl_std_median(:,2))];
%     ERPampl_skew = [ERPampl_skew, nanzscore(ERPampl_skew(:,2))];
%     ERPampl_kurt = [ERPampl_kurt, nanzscore(ERPampl_kurt(:,2))];
%     ERPampl_chi2stat = [ERPampl_chi2stat, nanzscore(ERPampl_chi2stat(:,2))];
% 
%     RR_ratio_std_mean_all_IC = [RR_ratio_std_mean_all_IC, nanzscore(RR_ratio_std_mean_all_IC(:,2))];
%     RR_skew_all = [RR_skew_all, nanzscore(RR_skew_all(:,2))];
%     RR_kurt_all = [RR_kurt_all, nanzscore(RR_kurt_all(:,2))];
%     RR_chi2stat_all_IC = [RR_chi2stat_all_IC, nanzscore(RR_chi2stat_all_IC(:,2))];
%     Rampl_ratio_std_mean_all_IC = [Rampl_ratio_std_mean_all_IC, nanzscore(Rampl_ratio_std_mean_all_IC(:,2))];
%     Rampl_skew_all_IC = [Rampl_skew_all_IC, nanzscore(Rampl_skew_all_IC(:,2))];
%     Rampl_kurt_all_IC = [Rampl_kurt_all_IC, nanzscore(Rampl_kurt_all_IC(:,2))];
%     Rampl_chi2stat_all_IC = [Rampl_chi2stat_all_IC, nanzscore(Rampl_chi2stat_all_IC(:,2))];
%     Qampl_ratio_std_mean_all_IC = [Qampl_ratio_std_mean_all_IC, nanzscore(Qampl_ratio_std_mean_all_IC(:,2))];
%     Qampl_skew_all_IC = [Qampl_skew_all_IC, nanzscore(Qampl_skew_all_IC(:,2))];
%     Qampl_kurt_all_IC = [Qampl_kurt_all_IC, nanzscore(Qampl_kurt_all_IC(:,2))];
%     Qampl_chi2stat_all_IC = [Qampl_chi2stat_all_IC, nanzscore(Qampl_chi2stat_all_IC(:,2))];
%     Sampl_ratio_std_mean_all_IC = [Sampl_ratio_std_mean_all_IC, nanzscore(Sampl_ratio_std_mean_all_IC(:,2))];
%     Sampl_skew_all_IC = [Sampl_skew_all_IC, nanzscore(Sampl_skew_all_IC(:,2))];
%     Sampl_kurt_all_IC = [Sampl_kurt_all_IC, nanzscore(Sampl_kurt_all_IC(:,2))];
%     Sampl_chi2stat_all_IC = [Sampl_chi2stat_all_IC, nanzscore(Sampl_chi2stat_all_IC(:,2))];
% 
%     % 2) Analyze distribution of RR_intervals and amplitude of R, S and Q peaks
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Determine the percentile we need to select the nb of IC chosen by user
%     nbr_IC_total = length(comp.label) - length(IC_bug) - length(IC_with_less_than_10_R_peak) - length(IC_to_not_analyze);
%     perc = nb_IC_wanted * 100 /nbr_IC_total;
% 
% 
%     % ERPampl
%     idx_score_ERPampl_median = find(ERPampl_median(:,2) < prctile(ERPampl_median(:,2), perc)  ); 
%     IC_ERPampl_median = ERPampl_median(idx_score_ERPampl_median,1);   
% 
%     idx_score_ERPampl_std = find(ERPampl_std(:,2) < prctile(ERPampl_std(:,2), perc)  ); 
%     IC_ERPampl_std = ERPampl_std(idx_score_ERPampl_std,1);   
% 
%     idx_score_ERPampl_skew = find(ERPampl_skew(:,2) < prctile(ERPampl_skew(:,2), perc)  ); 
%     IC_ERPampl_skew = ERPampl_skew(idx_score_ERPampl_skew,1);   
% 
%     idx_score_ERPampl_kurt = find(ERPampl_kurt(:,2) < prctile(ERPampl_kurt(:,2), perc)  ); 
%     IC_ERPampl_kurt = ERPampl_kurt(idx_score_ERPampl_kurt,1);   
% 
%     idx_score_ERPampl_std_median = find(ERPampl_std_median(:,2) < prctile(ERPampl_std_median(:,2), perc)  ); 
%     IC_ERPampl_std_median = ERPampl_std_median(idx_score_ERPampl_std_median,1);
% 
%     idx_score_ERPampl_chi2stat = find(ERPampl_chi2stat(:,2) < prctile(ERPampl_chi2stat(:,2), perc)  ); 
%     IC_ERPampl_chi2stat = ERPampl_chi2stat(idx_score_ERPampl_chi2stat,1);   
% 
% 
%     % RR_intervals
%     idx_score_std_mean = find(RR_ratio_std_mean_all_IC(:,2) < prctile(RR_ratio_std_mean_all_IC(:,2), perc)  );
%     IC_RR_std_mean = RR_ratio_std_mean_all_IC(idx_score_std_mean,1);
%     % figure; histogram(RR_ratio_std_mean_all_IC(:,2), 50);xlabel('RR std/mean'); ylabel('Nb of IC')
% 
%     idx_score_skew = find(RR_skew_all(:,2) < prctile(RR_skew_all(:,2), perc)  );
%     IC_RR_skew = RR_skew_all(idx_score_skew,1);
%     % figure; histogram(RR_skew_all(:,2), 20);xlabel('RR skewness'); ylabel('Nb of IC')
% 
%     idx_score_kurt = find(RR_kurt_all(:,2) < prctile(RR_kurt_all(:,2), perc)  );
%     IC_RR_kurt = RR_kurt_all(idx_score_kurt,1);
%     % figure; histogram(RR_kurt_all(:,2), 150);xlabel('RR kurtosis'); ylabel('Nb of IC')
% 
%     idx_chi2stat = find(RR_chi2stat_all_IC(:,2) < prctile(RR_chi2stat_all_IC(:,2), perc)  );
%     IC_RR_chi2stat = RR_chi2stat_all_IC(idx_chi2stat,1);
%     % figure; histogram(RR_chi2stat_all_IC(:,2), 80);xlabel('RR stat unif'); ylabel('Nb of IC')
% 
% 
%     % Amplitude of R peaks
%     idx_score_RR = find(Rampl_ratio_std_mean_all_IC(:,2) < prctile(Rampl_ratio_std_mean_all_IC(:,2), perc)  );
%     IC_Rampl_std_mean = Rampl_ratio_std_mean_all_IC(idx_score_RR,1);
%     % figure; histogram(Rampl_ratio_std_mean_all_IC(:,2), 150);xlabel('Rampl std/mean'); ylabel('Nb of IC')
% 
%     idx_skew = find(Rampl_skew_all_IC(:,2) < prctile(Rampl_skew_all_IC(:,2), perc)  );
%     IC_Rampl_skew = Rampl_skew_all_IC(idx_skew,1);
%     % figure; histogram(Rampl_skew_all_IC(:,2), 20);xlabel('Rampl skewness'); ylabel('Nb of IC')
% 
%     idx_score_kurt = find(Rampl_kurt_all_IC(:,2) < prctile(Rampl_kurt_all_IC(:,2), perc)  );
%     IC_Rampl_kurt = Rampl_kurt_all_IC(idx_score_kurt,1);
% 
%     idx_chi2stat = find(Rampl_chi2stat_all_IC(:,2) < prctile(Rampl_chi2stat_all_IC(:,2), perc)  );
%     IC_Rampl_chi2stat = Rampl_chi2stat_all_IC(idx_chi2stat,1);
%     % figure; histogram(Rampl_chi2stat_all_IC(:,2), 100);xlabel('Rampl stat unif'); ylabel('Nb of IC')
% 
% 
%     % Amplitude of Q peaks
%     idx_score_Q = find(Qampl_ratio_std_mean_all_IC(:,2) < prctile(Qampl_ratio_std_mean_all_IC(:,2), perc)  ); 
%     IC_Qampl_std_mean = Qampl_ratio_std_mean_all_IC(idx_score_Q,1);
% 
%     idx_skew = find(Qampl_skew_all_IC(:,2) < prctile(Qampl_skew_all_IC(:,2), perc)  );
%     IC_Qampl_skew = Qampl_skew_all_IC(idx_skew,1);
% 
%     idx_score_kurt = find(Qampl_kurt_all_IC(:,2) < prctile(Qampl_kurt_all_IC(:,2), perc)  );
%     IC_Qampl_kurt = Qampl_kurt_all_IC(idx_score_kurt,1);
% 
%     idx_chi2stat = find(Qampl_chi2stat_all_IC(:,2) < prctile(Qampl_chi2stat_all_IC(:,2), perc)  );
%     IC_Qampl_chi2stat = Qampl_chi2stat_all_IC(idx_chi2stat,1);
% 
% 
% 
%     % Amplitude of S peaks
%     idx_score_S = find(Sampl_ratio_std_mean_all_IC(:,2) < prctile(Sampl_ratio_std_mean_all_IC(:,2), perc)  );  
%     IC_Sampl_std_mean = Sampl_ratio_std_mean_all_IC(idx_score_S,1);
% 
%     idx_skew = find(Sampl_skew_all_IC(:,2) < prctile(Sampl_skew_all_IC(:,2), perc)  );
%     IC_Sampl_skew = Sampl_skew_all_IC(idx_skew,1);
% 
%     idx_score_kurt = find(Sampl_kurt_all_IC(:,2) < prctile(Sampl_kurt_all_IC(:,2), perc)  );
%     IC_Sampl_kurt = Sampl_kurt_all_IC(idx_score_kurt,1);
% 
%     idx_chi2stat = find(Sampl_chi2stat_all_IC(:,2) < prctile(Sampl_chi2stat_all_IC(:,2), perc)  );
%     IC_Sampl_chi2stat = Sampl_chi2stat_all_IC(idx_chi2stat,1);
% 
% 
% 
%     % 3) Find IC that meet bpm and SignalAmpl range conditions
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % Bpm
%     idx_bpm = intersect(  find((bpm_all_IC(:,2) > bpm_min)),  find(bpm_all_IC(:,2) < bpm_max) );
%     IC_bpm_ok = bpm_all_IC(idx_bpm,1);
% 
%     % SignalAmpl range
%     idx_SignalAmpl_range = find(SignalAmpl_range_all_IC(:,2) < threshold_regularity_signal_minmax);
%     IC_SignalAmpl_range_ok = SignalAmpl_range_all_IC(idx_SignalAmpl_range,1);
% 
% 
% 
%     %% Find heart IC (METHOD A: WITH ZSCORE)
% 
%     % 1) Gather all zscores in a table
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Zscore_ALL = [];
%     for IC_iter = 1:size(RR_ratio_std_mean_all_IC,1)
%         IC_name = RR_ratio_std_mean_all_IC(IC_iter, 1);
%         Zscore_ALL = [Zscore_ALL; ...
%             IC_name, ERPampl_median(IC_iter,3), ERPampl_std(IC_iter,3), ERPampl_skew(IC_iter,3), ERPampl_kurt(IC_iter,3), ERPampl_chi2stat(IC_iter,3), ERPampl_std_median(IC_iter,3),...
%             RR_ratio_std_mean_all_IC(IC_iter,3), RR_skew_all(IC_iter,3), RR_kurt_all(IC_iter,3), RR_chi2stat_all_IC(IC_iter,3), ...
%             Rampl_ratio_std_mean_all_IC(IC_iter,3), Rampl_skew_all_IC(IC_iter,3), Rampl_kurt_all_IC(IC_iter,3), Rampl_chi2stat_all_IC(IC_iter,3), ...
%             Qampl_ratio_std_mean_all_IC(IC_iter,3), Qampl_skew_all_IC(IC_iter,3), Qampl_kurt_all_IC(IC_iter,3), Qampl_chi2stat_all_IC(IC_iter,3), ...
%             Sampl_ratio_std_mean_all_IC(IC_iter,3), Sampl_skew_all_IC(IC_iter,3), Sampl_kurt_all_IC(IC_iter,3), Sampl_chi2stat_all_IC(IC_iter,3)];
%     end
%     Zscore_ALL = array2table(Zscore_ALL);
%     Zscore_ALL.Properties.VariableNames = {'IC', 'z_ERPampl_median', 'z_ERPampl_std', 'z_ERPampl_skew', 'z_ERPampl_kurt', 'z_ERPampl_chi2stat', 'z_ERPampl_std_median', ...
%             'z_RR_ratio_std_mean', 'z_RR_skew', 'z_RR_kurt', 'z_RR_chi2stat', ...
%             'z_Rampl_ratio_std_mean', 'z_Rampl_skew', 'z_Rampl_kurt', 'z_Rampl_chi2stat', ...
%             'z_Qampl_ratio_std_mean', 'z_Qampl_skew', 'z_Qampl_kurt', 'z_Qampl_chi2stat', ...
%             'z_Sampl_ratio_std_mean', 'z_Sampl_skew', 'z_Sampl_kurt', 'z_Sampl_chi2stat'};
% 
%     % 2) Define coeff of zscore (as for some metrics the heart IC will have the lowest zscore, but will have the highest zscore for other metrics)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % var_coeff = [0, -1, repmat(1,1,8), -1,-1,1,1,  -1,-1,1,1];
%     % var_coeff = array2table(var_coeff);
%     % var_coeff.Properties.VariableNames = {'NA', 'coeff_ERP_ampl', 'coeff_RR_ratio_std_mean', 'coeff_RR_skew', 'coeff_RR_kurt', 'coeff_RR_chi2stat', ...
%     %         'coeff_Rampl_ratio_std_mean', 'coeff_Rampl_skew', 'coeff_Rampl_kurt', 'coeff_Rampl_chi2stat', ...
%     %         'coeff_Qampl_ratio_std_mean', 'coeff_Qampl_skew', 'coeff_Qampl_kurt', 'coeff_Qampl_chi2stat', ...
%     %         'coeff_Sampl_ratio_std_mean', 'coeff_Sampl_skew', 'coeff_Sampl_kurt', 'coeff_Sampl_chi2stat'};
% 
% 
%     % 3) Compute ponderated zscore for each IC (i.e., sum of (zscore{i} x coeff{i}))
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Z_score_sum_all = [];
%     col_of_interest = [2:4,8:9]; % (2:6 = ERPampl_median, std, skew, kurt, chi2stat)
%     for IC_iter = 1:size(RR_ratio_std_mean_all_IC,1)
%         Z_score_sum = 0;
% 
%         for col_iter = 1:length(col_of_interest)
%             col_iter = col_of_interest(col_iter);
%             Z_score_sum = Z_score_sum + table2array(Zscore_ALL(IC_iter,col_iter)); % * table2array(var_coeff(1,col_iter));
%         end
% 
%         % % Add the best zscores of Q or S peaks
%         %     % Compute Q and S features (Q: col 15 to 18, S: col 19 to 22)
%         %     Q_iter = 0;
%         %     S_iter = 0;
%         %     for col_iter = 15:18
%         %         Q_iter = Q_iter + table2array(Zscore_ALL(IC_iter,col_iter)); % * table2array(var_coeff(1,col_iter));
%         %         S_iter = S_iter + table2array(Zscore_ALL(IC_iter,col_iter+4)); % * table2array(var_coeff(1,col_iter+4));
%         %     end
%         % 
%         %     % Add max zscore (Q or S features)
%         %     Z_score_sum = Z_score_sum + min(Q_iter, S_iter);
% 
%         % Gather with the Z_score_sum of the other IC
%         Z_score_sum_all = [Z_score_sum_all; Z_score_sum];
%     end
% 
%     Z_score_sum_all = array2table(Z_score_sum_all);
%     Z_score_sum_all.Properties.VariableNames = {'Ponderated_sum_zscore'};
% 
%     Zscore_ALL = [Zscore_ALL, Z_score_sum_all];
% 
%     % 5) Compute threshold to identify heart IC
%     thresh_zscore = nanmean(Zscore_ALL.Ponderated_sum_zscore) - threshold_std_method2*nanstd(Zscore_ALL.Ponderated_sum_zscore);
%     min(Zscore_ALL.Ponderated_sum_zscore)
%     find(Zscore_ALL.Ponderated_sum_zscore < thresh_zscore);
% 
%     % 6) Contribution of each metric
%         % Find the IC with lowest z-score_sum
%     [min_zscore, potential_IC] = min(Zscore_ALL.Ponderated_sum_zscore);
%         % Extract all its z-scores
%     z_score_potential_heart_IC = Zscore_ALL(potential_IC,2:end-1);
%         % Sort its z-scores from lowest to highest
%     [val, idx_sorting] = sort(table2array(z_score_potential_heart_IC));
%     z_score_potential_heart_IC = z_score_potential_heart_IC(:, idx_sorting); 
% 
%         % Extract the strength of the top 3 metrics
%     percentage_of_zscore_with_only_top3_metrics = 100 * sum(table2array(z_score_potential_heart_IC(1,1:3))) / Zscore_ALL.Ponderated_sum_zscore(potential_IC); % and not "sum(table2array(z_score_potential_heart_IC(1,1:end)))" bc we kept only the best zscore of either Q or S
% 
% 
% 
%     % LAST) Remove some IC based on bpm and SignalAmp
%     Zscore_ALL.Ponderated_sum_zscore(find(~ismember(1:size(Zscore_ALL,1),IC_bpm_ok))) = NaN;
%     Zscore_ALL.Ponderated_sum_zscore(find(~ismember(1:size(Zscore_ALL,1),IC_SignalAmpl_range_ok))) = NaN;
% 
% 
%     %% Find heart IC (METHOD B: With BINARY score)
% 
%     % 1) Method 1: Compute final BINARY score and find IC with score above threshold
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %%%%%%%%%%%%%%%%%%%%%%% USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % // To modifiy if you change the nb of conditions below \\
%     nbr_condition = 5;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     IC_heart_method1 = [];
%     score_other_IC = [];
%     all_score = [];
% 
%     table_IC_cond = [];
% 
%     for comp_iter = 1:length(comp.label)
%         % Case IC_bug
%         if ismember(comp_iter, IC_bug)
%             score_other_IC = [score_other_IC, NaN];
%             all_score = [all_score; comp_iter, NaN];
%             table_IC_cond = [table_IC_cond; comp_iter, repmat(NaN,1,nbr_condition+2)]; % +2 for bpm and IC_SignalAmpl_range_ok
%         else
% 
%             % Check if the IC iter meet the different conditions
%             %ERP_ampl (QRS)
%             cond01 = ismember(comp_iter, IC_ERPampl_median);
%             cond02 = ismember(comp_iter, IC_ERPampl_skew);
%             % cond03 = ismember(comp_iter, IC_ERPampl_kurt);
%             % cond04 = ismember(comp_iter, IC_ERPampl_chi2stat);
%             % cond05 = ismember(comp_iter, IC_ERPampl_std);
%             cond06 = ismember(comp_iter, IC_ERPampl_std_median);
% 
% 
% 
% 
%             % RR interval conditions
%             cond1 = ismember(comp_iter, IC_RR_std_mean);
%             cond2 = ismember(comp_iter, IC_RR_skew);
%             % cond3 = ismember(comp_iter, IC_RR_kurt);
%             % cond4 = ismember(comp_iter, IC_RR_chi2stat);
%             % Rampl conditions
%             % cond5 = ismember(comp_iter, IC_Rampl_std_mean);
%             % cond6 = ismember(comp_iter, IC_Rampl_skew);
%             % cond7 = ismember(comp_iter, IC_Rampl_kurt);
%             % cond8 = ismember(comp_iter, IC_Rampl_chi2stat);
%             % % Q and S ampl conditions
%             % cond9 = ismember(comp_iter, IC_Qampl_std_mean) || ismember(comp_iter, IC_Sampl_std_mean);
%             % cond10 = ismember(comp_iter, IC_Qampl_skew) || ismember(comp_iter, IC_Sampl_skew);
%             % cond11 = ismember(comp_iter, IC_Qampl_kurt) || ismember(comp_iter, IC_Sampl_kurt);
%             % cond12 = ismember(comp_iter, IC_Qampl_chi2stat) || ismember(comp_iter, IC_Sampl_chi2stat);
% 
%             % Proportion of conditions met
%             % prop_cond_ok = (cond01 + cond1 + cond2 + cond3 + cond4 + cond5 + cond6 + cond7 + cond8 + cond9+cond10+cond11+cond12 ) / nbr_condition ;
%             prop_cond_ok = (cond01+cond02+cond06 + cond1+cond2) / nbr_condition ;
% 
%             all_score = [all_score; comp_iter, prop_cond_ok];
% 
%             % Specific conditions that heart IC MUST meet (bpm and SignalAmplRange)
%             cond_bpm = ismember(comp_iter, IC_bpm_ok);
%             cond_SignalAmpl_range = ismember(comp_iter, IC_SignalAmpl_range_ok);
% 
%             % Gather info with other IC
%             % table_IC_cond = [table_IC_cond; comp_iter, cond01, cond1, cond2, cond3, cond4, cond5, cond6, cond7, cond8, cond9, cond10, cond11, cond12, cond_bpm, cond_SignalAmpl_range]; 
%             table_IC_cond = [table_IC_cond; comp_iter, cond01,cond02,cond06, cond1,cond2, cond_bpm, cond_SignalAmpl_range]; 
% 
%             % Determine if it's a heart IC
%             if (prop_cond_ok >= threshold_cond_IC_method1) && (cond_bpm == 1) && (cond_SignalAmpl_range == 1)
%                 IC_heart_method1 = [IC_heart_method1; comp_iter, prop_cond_ok];
%                 score_other_IC = [score_other_IC, NaN];
%             else
%                 score_other_IC = [score_other_IC, prop_cond_ok];
%             end
%         end
%     end
% 
%     % Give name to table_IC_cond col
%     table_IC_cond = array2table(table_IC_cond);
%     % table_IC_cond.Properties.VariableNames = {'IC', 'ERP_ampl', 'RR_std_mean', 'RR_skew', 'RR_kurt', 'RR_chi2stat', ...
%     %     'Rampl_std_mean', 'Rampl_skew', 'IC_Rampl_kurt', 'IC_Rampl_chi2stat', ...
%     %     'Q_or_Sampl_std_mean', 'Q_or_Sampl_skew', 'Q_or_Sampl_kurt', 'Q_or_Sampl_chi2stat',...
%     %     'bpm_physiologic', 'SignalAmpl_range_ok'};
%     table_IC_cond.Properties.VariableNames = {'IC', 'ERP_ampl_median', 'ERP_ampl_skew', 'ERP_ampl_std_median', ...
%         'RR_std_mean', 'RR_skew', ...
%         'bpm_physiologic', 'SignalAmpl_range_ok'};
% 
%     % Extract IC_num, bpm and score for the heart IC
%     if isempty(IC_heart_method1) == 0
%         IC_heart_method1_num = IC_heart_method1(:,1);
%         IC_heart_method1_score = IC_heart_method1(:,2);
%         IC_heart_method1_bpm = bpm_all_IC(ismember(bpm_all_IC(:,1), IC_heart_method1_num), 2);
%     else
%         IC_heart_method1_num = [];
%         IC_heart_method1_score = [];
%         IC_heart_method1_bpm = [];
%     end
% 
%     % To have an idea of the proportion of conditions met by non-heart IC
%     Max_score_other_IC = max(score_other_IC);
% 
%     % For sanity check
%     Nbr_sec_ICA = size(comp.trial,2);
% 
% 
%     % 2) Method 2: find heart IC based on the distribution of scores
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Candidates = IC with score > mean + n*std
%     IC_heart_method2_candidates = all_score(find( all_score(:,2) > nanmean(all_score(:,2)) + threshold_std_method2*nanstd(all_score(:,2)) ) ,1);
% 
%     % For each candidate IC, check if 1) bpm physiological and 2) the SignalAmplRange
%     IC_heart_method2 = [];
%     for i = 1:length(IC_heart_method2_candidates)
%         iter = IC_heart_method2_candidates(i);
%         if ismember(iter, IC_bpm_ok) && ismember(iter, IC_SignalAmpl_range_ok) && ismember(iter,IC_heart_method2_candidates)
%             IC_heart_method2 = [IC_heart_method2; iter];
%         end
%     end
% 
% 
%     %% Summarize info in a table
% 
%     % Create a table for recording iter
%     % OLD VERSION (without IC n° for bpm_all_IC and SignalAmpl_range_all_IC)
%     % table_cardiac_IC = table({IC_bug}, {bpm_all_IC(:,2)}, {IC_bpm_ok}, {SignalAmpl_range_all_IC(:,2)}, {IC_SignalAmpl_range_ok},...
%     %     {IC_RR_std_mean}, {IC_RR_skew}, {IC_RR_kurt}, {IC_RR_chi2stat},  ...
%     %     {IC_Rampl_std_mean}, {IC_Rampl_skew}, {IC_Rampl_kurt}, {IC_Rampl_chi2stat}, ...
%     %     {IC_heart_method2}, {IC_heart_method1_num}, {IC_heart_method1_score}, {all_score}, {score_other_IC}, {Max_score_other_IC}, {IC_heart_method1_bpm}, {Nbr_sec_ICA});
% 
%     table_cardiac_IC = table({IC_bug}, {bpm_all_IC}, {IC_bpm_ok}, {SignalAmpl_range_all_IC}, {IC_SignalAmpl_range_ok},...
%         {IC_ERPampl_median}, {IC_ERPampl_std},{IC_ERPampl_skew},{IC_ERPampl_kurt}, {IC_ERPampl_chi2stat}, {IC_ERPampl_std_median},...
%         {IC_RR_std_mean}, {IC_RR_skew}, {IC_RR_kurt}, {IC_RR_chi2stat},  ...
%         {IC_Rampl_std_mean}, {IC_Rampl_skew}, {IC_Rampl_kurt}, {IC_Rampl_chi2stat}, ...
%         {IC_heart_method2}, {IC_heart_method1_num}, {IC_heart_method1_score}, {all_score}, {score_other_IC}, {Max_score_other_IC}, {IC_heart_method1_bpm}, {Nbr_sec_ICA});
% 
% 
%     % Rename table columns
%     table_cardiac_IC.Properties.VariableNames = col_names;
% 
% 
% 
%     %% Extract the cardiac IC
% 
%     % Extract the cardiac IC depending on the chosing method
%     if strcmp(method_chosen, 'absolute_threshold')
%         heart_IC = table_cardiac_IC.IC_heart_method1{1}';
%     elseif strcmp(method_chosen, 'mean_std')
%         heart_IC = table_cardiac_IC.IC_heart_method2{1}';
%     else
%         error("Error: method_chosen must take the value 'absolute_threshold' or 'mean_std'")
%     end
% 
% 
%     % Extract its z-scores
%     z_score_heart_IC = Zscore_ALL(heart_IC,2:end-1);
%     % Sort its z-scores from lowest to highest
%     [val, idx_sorting] = sort(table2array(z_score_heart_IC));
%     z_score_heart_IC_sort = z_score_heart_IC(:, idx_sorting);
% 
% 
%     %% Plot and save time course of heart IC
%     if plot_heart_IC == 1
% 
%         % 1) Plot DISTRIBUTION PROP CONDITION OK
%         % threshold_heart1 = threshold_cond_IC_method1;
%         % threshold_heart2 = nanmean(all_score(:,2)) + threshold_std_method2*nanstd(all_score(:,2));
%         % figure
%         % plot(all_score(:,2),'.','linestyle','none', 'MarkerSize', 15);
%         % hold on
%         % xlim([0 length(all_score(:,2))+1]);
%         % xl = xlim;
%         % % yl = ylim;
%         % set(gca,'ylimMode','manual');
%         % plot(xl, [threshold_heart1, threshold_heart1], 'r--', 'LineWidth', 1);
%         % plot(xl, [threshold_heart2, threshold_heart2], 'b--', 'LineWidth', 1);
%         % % Emphasize the rejected IC
%         % hold on
%         % emphasizeX = [heart_IC];    % X-coordinates of points to emphasize
%         % emphasizeY = all_score(heart_IC,2);  % Corresponding Y-coordinates
%         % plot(emphasizeX, emphasizeY, '.','linestyle','none', 'MarkerSize', 25, 'Color', 'r');
%         % title(['Prop conditions ok (red: absolute thresh, blue:mean score*' num2str(threshold_std_method2) 'std )'])
%         % ylabel('Heart score (using pop_iclabel)');
%         % xlabel('Components');
%         % % %%%%% SAVE fig as .png %%%%%%%%
%         % filename = strcat(cfg.path_output, '/Distribution_score/', file_info, '_ECG_IC_label_fct_find_cardiac.png');
%         % saveas(gcf, filename);
%         % close(gcf);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % 2) Plot all IC time course (cardiac and non-cardiac in separate folders)
%         % for comp_iter = 1:size(comp.label,1)
%         % 
%         %     if ismember(comp_iter, heart_IC)
%         %         folder_final = 'heart';
%         %     else
%         %         folder_final = 'non_heart';
%         %     end
%         % 
%         %     figure
%         %     IC_to_plot = comp_iter; %heart_IC(comp_iter);
%         %     plot([1:size(IC_timecourse,2)], IC_timecourse(IC_to_plot,:))
%         %     title(['Time course (IC ' num2str(IC_to_plot) ')'])
%         %     % %%%%% SAVE fig as .png %%%%%%%%
%         %     filename = strcat(cfg.path_output, '/IC_timecourse/', folder_final, '/', file_info, '_ECG_fct_find_cardiac_time_course_heart_IC ', num2str(IC_to_plot), '.png');
%         %     saveas(gcf, filename);
%         %     close(gcf);
%         %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 
%         % end
% 
% 
% 
            figure(94480); clf;
            set(gcf,'WindowState', 'maximized'); % Open a maximized figure window

%         %%%%%%%%%%%%%%% Several plots per window %%%%%%%%%%%%%
        nbr_plot_per_window = 9; % If it is changed, change also subplot parameters
        nbr_IC = length(comp.label);

        for i = 1 : nbr_plot_per_window : nbr_IC


            IC_start = i;
            IC_stop =  min(i + nbr_plot_per_window - 1, nbr_IC); % Make sure we don't go beyond n
            
            figure(94480); clf;
            for IC_to_plot = IC_start:IC_stop
                position_idx = mod(IC_to_plot, nbr_plot_per_window); % 10 --> 1, 11 --> 2...
                if position_idx == 0
                    position_idx = 9;
                end
                subplot(3, 3, position_idx); 

                % Generate a plot
                timetoplot = [0 20]; % <<- parametrize?
                t = timepts(timetoplot, comp.time{1});
                plot(comp.time{1}(t), comp.trial{1}(IC_to_plot,t), ifelse(ismember(IC_to_plot, heart_IC),'r',''));

                % Add a title and grid
                title(['IC ' num2str(IC_to_plot)])
                grid on;

                % Optional: Label axes (only if necessary)
                xlabel('Time (s)');
                ylabel('Voltage');
            end
            % %%%%% SAVE fig as .png %%%%%%%%
            filename = strcat(cfg.path_output, '/IC_timecourse/', file_info, '_IC ', num2str(IC_start), '_to_IC_', num2str(IC_stop), '.png');
            mymkdir(fileparts(filename))
            saveas(gcf, filename);
            close(gcf);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
% 
% 
% 
%     end
% 
% 
% end
% 
% 
% %% NEW FOR TESTING
% 
% % Distribution of ERP_ampl
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % ERPampl_median
% % hist(ERPampl_median(:,2))
% % nanmean(ERPampl_median(:,2)) + 2*nanstd(ERPampl_median(:,2))
% 
% 
% % Plot ERP
% %%%%%%%%%%%
% % IC_to_plot = 2; 
% % figure; 
% % for i = 1:size(ERP_details_for_each_IC{1, IC_to_plot}, 1)
% %     hold on
% %     plot(ERP_details_for_each_IC{1,IC_to_plot}(i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
% % end
% % plot(ERP_average_for_each_IC(IC_to_plot,:), 'Color', 'red', 'LineWidth', 2)
% 
% 
% % Plot ERP of some IC
% %%%%%%%%%%%%%%%%%%%%%
% % for comp_iter = 1:4
% %     figure
% %     plot([1:length(ERP_average_for_each_IC)], ERP_average_for_each_IC(comp_iter,:))
% %     title(num2str(comp_iter))
% % end
% 
% 
% 
% 
% %% Correlation matrix
% 
% min_scale = -1;
% max_scale = 1;
% 
% % Gather all features of all metric
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL_ALL_FEATURES = [ERPampl_median(:,2), ERPampl_std(:,2),  ERPampl_skew(:,2),  ERPampl_kurt(:,2),  ERPampl_std_median(:,2), ...
%     RR_ratio_std_mean_all_IC(:,2), RR_skew_all(:,2), RR_kurt_all(:,2), RR_chi2stat_all_IC(:,2), ...
%     Rampl_ratio_std_mean_all_IC(:,2), Rampl_skew_all_IC(:,2), Rampl_kurt_all_IC(:,2), Rampl_chi2stat_all_IC(:,2), ...
%     Qampl_ratio_std_mean_all_IC(:,2), Qampl_skew_all_IC(:,2), Qampl_kurt_all_IC(:,2), Qampl_chi2stat_all_IC(:,2), ...
%     Sampl_ratio_std_mean_all_IC(:,2), Sampl_skew_all_IC(:,2), Sampl_kurt_all_IC(:,2), Sampl_chi2stat_all_IC(:,2)];
% 
% 
% % Compute corr matrix
% %%%%%%%%%%%%%%%%%%%%%
% correlationMatrix = corr(FINAL_ALL_FEATURES, 'Rows', 'complete');
% 
% 
% % List variable names for corr matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varNames = {'ERPampl_median', 'ERPampl_std', 'ERPampl_skew', 'ERPampl_kurt', 'ERPampl_std_median',...
%     'RR_ratio_std_mean', 'RR_skew', 'RR_kurt', 'RR_chi2stat', ...
%     'Rampl_ratio_std_mean', 'Rampl_skew', 'Rampl_kurt', 'Rampl_chi2stat', ...
%     'Qampl_ratio_std_mean', 'Qampl_skew', 'Qampl_kurt', 'Qampl_chi2stat', ...
%     'Sampl_ratio_std_mean', 'Sampl_skew', 'Sampl_kurt', 'Sampl_chi2stat'};
% 
% 
% % Plot corr matrix
% %%%%%%%%%%%%%%%%%%
% % figure
% % imagesc(correlationMatrix);
% % colorbar;
% % colormap(jet);
% % xticks(1:length(varNames)); xticklabels(varNames);
% % yticks(1:length(varNames)); yticklabels(varNames);
% % % Set color limits (e.g., [-1, 1] for correlation)
% % caxis([min_scale max_scale]);                          % Define color bar limits
% % colorbar;                               % Add color bar
% % colormap(jet);                          % Use a jet colormap
% % title(['Correlation Matrix (' strrep(cfg.file_info, '_', ' ') ')']);
% 
% 
% % Annotate each cell with its r-value
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [nRows, nCols] = size(correlationMatrix); % Get the matrix dimensions
% % for i = 1:nRows
% %     for j = 1:nCols
% %         text(j, i, num2str(correlationMatrix(i, j), '%.2f'), ... % Format values to 2 decimal places
% %             'HorizontalAlignment', 'center', ...
% %             'Color', 'white', ...         % Text color
% %             'FontSize', 10);             % Adjust font size
% %     end
% % end
% 
% % Save
% %%%%%%
% % filename = strcat([cfg.path_output, '/corr_matrices/Corr_matrix_' cfg.file_info '.png']);
% % saveas(gcf, filename);
% % close(gcf);
% 
% 
% %% Extract all features for ROC
% 
% FINAL_ALL_FEATURES_table = array2table(FINAL_ALL_FEATURES);
% FINAL_ALL_FEATURES_table.Properties.VariableNames = varNames;
% 
% 
% 
% %% Reduce the number of outputs
% 
% output_for_zscore_corMatrix_ROC = [];
% output_for_zscore_corMatrix_ROC.z_score_heart_IC = z_score_heart_IC;
% output_for_zscore_corMatrix_ROC.FINAL_ALL_FEATURES_table = FINAL_ALL_FEATURES_table;
% output_for_zscore_corMatrix_ROC.correlationMatrix = correlationMatrix;
% 
% 
% 
% output_for_user = [];
% output_for_user.heart_IC = heart_IC;
% output_for_user.table_IC_cond = table_IC_cond;
% output_for_user.table_cardiac_IC = table_cardiac_IC;

end


