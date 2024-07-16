function [ecg_f] = preprocess_window_ecg(ecg,fs) %#codegen
%% WHAT DOES THIS FUNCTION DO?
% THIS FUNCTION FILTERS AND PRE-PROCESSES THE UNFILTERED RAW ECG SIGNAL,

% NORMALLY WE APPLY A LOWPASS FILTER AND A NOTCH FILTER.

% AS FOR THE WINDOWING CODE, THE CENTRE R PEAK IS COMPUTED. THE WONDOW
% SIZE, IS 5 TIMES THE MEAN RR INTERVAL.

% THE START INDEX IS FROM R(CENTRE) - W/2 -> ROUND OFF TO NEAREST INTEGER
% VALUE.

% WINDOWING IS DONE, TO ENSURE, THAT THE ECG WAVE IS SYMMETRIC AND THEN THE
% P AND T PEAKS CAN BE MARKED AND DIFFERENTIATED.
%% AUTHOR- 
% COMPILED AND MAINTAINED BY-
% ROHAN SANGHAVI.
%% FILTERING

b = [-0.000187830184206340,-0.000205456041853291,-8.64057295452613e-05,0.000234938106149921,0.000787322828783533,0.00153305158366326,0.00234705727203288,0.00301420392436277,0.00325208908216934,0.00276166534082112,0.00130116641212287,-0.00122845080648299,-0.00470423015041083,-0.00872851920245046,-0.0126215652355866,-0.0154716371701053,-0.0162444264117259,-0.0139408664255059,-0.00778038794719607,0.00262217964633411,0.0171209674455062,0.0349757093199774,0.0548809558645369,0.0750922598049386,0.0936356537755463,0.108570111841377,0.118261039548750,0.121618807015993,0.118261039548750,0.108570111841377,0.0936356537755463,0.0750922598049386,0.0548809558645369,0.0349757093199774,0.0171209674455062,0.00262217964633411,-0.00778038794719607,-0.0139408664255059,-0.0162444264117259,-0.0154716371701053,-0.0126215652355866,-0.00872851920245046,-0.00470423015041083,-0.00122845080648299,0.00130116641212287,0.00276166534082112,0.00325208908216934,0.00301420392436277,0.00234705727203288,0.00153305158366326,0.000787322828783533,0.000234938106149921,-8.64057295452613e-05,-0.000205456041853291,-0.000187830184206340];
a = 1;
ecg1 = filter(b,a,ecg);
ecg_f = baseline_remove(ecg1);    
ecg_f = ecg_f./max(ecg_f);
ecg_f = ecg_f(100:end);


[~,locs_R,delay]=pan_tompkin2(ecg_f,fs); %more accurate.
search_R = round(0.08*fs); % 80ms windowing time;
locs_Rn = locs_R - round(delay);

check_loc2 = find(locs_Rn>=length(ecg_f) | locs_Rn - search_R <=0 | locs_Rn + search_R >=length(ecg_f));
locs_Rn((check_loc2)) = [];
locs_R1 = locs_Rn + search_R;
locs_R2 = locs_Rn - search_R;


locs_Rf = zeros(1,length(locs_Rn));
for k = 1:length(locs_Rn)
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %%%%%%%%%% OLD (bug if two timepoint have exactly the same value
    % locs_Rf(k) = find(ecg_f == max(ecg_f(locs_R2(k):locs_R1(k))));
    %%%%%%%%%% PIERRE
    % Find the candidates for the peak
    temp_PC =  find(ecg_f == max(ecg_f(locs_R2(k):locs_R1(k))));
    % If several values, find the one closer to the P peak
    if length(temp_PC) > 1
        % Compute dist to the P peak
        dist_with_peak = [];
        for iter_PC = 1:length(temp_PC)
            value_iter_PC = temp_PC(iter_PC);
            dist_with_peak = [dist_with_peak, abs(locs_Rf(k-1) - value_iter_PC)];
        end
        % Find the min dist locs_Pf
        [minValue_PC, minIndex_PC] = min(dist_with_peak);
        % Add it to locs_Pf
        locs_Tf(k) = temp_PC(minIndex_PC);
    else
        locs_Tf(k) =  temp_PC;
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
end
mean_RR = mean(diff(locs_Rf))/fs; 
%% WINDOWING
Number = change_odd(length(locs_Rf) - round(length(locs_Rf)/5)); % change on convenience. In the paper it is given 
% change_odd(length(locs_Rf) - 5); the reason being the signal is quite
% long. Here it is short. Please change accordingly when signal is long or
% short. This was finely tuned to the application which my group had
% designed.
% if one does not wish to keep tuning use what I have given here. Thank
% you!!!

wl = round(mean_RR*fs)*Number;
idx = locs_R(round(length(locs_R)/2)) - round(wl/2); %this is done to ensure that the ecg signal
%starts symetrically about the the R peak as often as possible
%so at centre R peak, - 1/2 the window length so symmetry is
%assured.(VVVVVIMP)

%%%%%%%%%%% PIERRE CHAMPETIER %%%%%%%%%%%
% For some IC, I have errors bc ::
% 1) idx < 0 or 2) wl + idx - 1 > size(ecg_f,2), or 3) idx = NaN
if idx < 0 ||  wl + idx - 1 > size(ecg_f,2)
    idx = 1;
end
% 
% if isnan(wl)
%     wl = size(ecg_f,2);
% end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ecg_f = ecg_f(bsxfun(@plus,idx,0:wl-1));  %extracted indexes %%% COMMENTED BY PIERRE
ecg_f = detrend(ecg_f,10);
ecg_f = ecg_f./max(ecg_f);
end 