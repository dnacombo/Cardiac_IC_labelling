function [locs_Pf,locs_Q,locs_Rf,locs_Sf,locs_Tf,varargout] = compute_fudicial_peaks_live17_c(ecg_f,fs,plot_data) %#codegen
%% WHAT DOES THIS FUNCTION DO?
% THIS FUNCTION COMPUTES IN REAL TIME THE P,Q,R,S AND T PEAKS OF AN ECG
% SIGNAL.
% R PEAKS- PAN TOMPKIN CODE IS USED DIRECTLY SERACHING IS TAKEN FROM NIMA
% AALIZADE (REF).

% S PEAKS - WE SEARCH AFTER THE R PEAKS TO FIND THE S PEAKS USING A SERACH
% OFFSET OF 0.1*fs.

% Q PEAKS- WE SEARCH BEFORE THE R PEAKS TO FIND THE Q PEAKS. SERAH OFFSET
% IS 0.08*FS (80ms)

% P AND T PEAKS - IN THE UPDATED VERSION THE THRESHOLDING IS COMPLETELY
% ELIMIMTATED USING THE IDEA AND CODE PROVIDED IN 
% WE SEARCH 1/2 THE MEAN RR DISTANCE BEFORE R PEAK AND AFTER S PEAK. THIS
% WAY THE SEARCH OFFSET KEEPS VARYING AND IT NEVER CROSSES THE R PEAK ALSO
% (ORIGINAL).

% THE PROPOSED WINDOWING CODE HELPS TO GET SYMMETRY IN THE ECG WAVEFORM.
% ALSO THERE IS A ADAPTIVE THRESHOLD TO COMPUTE P PEAKS AND AVOID FLUTTER.
%% MAIN REFERENCE
% https://in.mathworks.com/matlabcentral/fileexchange/66098-ecg-p-qrs-t-wave-detecting-matlab-code.
%% FINAL STEPS (ALGORITHM).

% STEP 1-> SEARCH FOR Q,S,P AND T PEAKS FROM THE R PEAK (ONE AFTER
% ANOTHER). ( REFERRED-> NIMA AALIZADE)

% STEP2 -> ADAPTIVE THRESHOLDING TO ELIMINATE FAKE P PEAKS MARKED IN AF &
% FLUTTER
% DONE!!!!!
%% EXTRA INFORMATION
% THE P WAVE IS ACCEPTED ONLY IF IT IS BETWEEN 1/10 AND 1/4 OF THE
% CORRESPONDING R WAVE.
% R AMPLITUDE -> 2.5mV(amplitude)
% P AMPLITUDE -> 0.25mV (min)
% so 10 ( 20 if not detected 2nd chance).
% p waves cannot be extremely high hence 1/4.
% IF IT IS ACCEPTED THE THRESHOLDS ARE FURTHER STRETCHED TO INCLUDE
% BASELINE WANDER.
% OR ELSE REJECTED.
%% AUTHOR

% COMPILED AND MAINTAINED BY-
% ROHAN SANGHAVI.
%% PARTIAL CREDIT.> NIMA AALIZADE.
% P AND T PEAK COMPUTATON REFERRED FROM-
% https://in.mathworks.com/matlabcentral/fileexchange/66098-ecg-p-qrs-t-wave-detecting-matlab-code
% https://in.mathworks.com/matlabcentral/fileexchange/66098-ecg-p-qrs-t-wave-detecting-matlab-code
%% PAN TOMPKIN CODE FOR COMPUTATION OF R PEAKS ( CODE REFERRED DIRECTLY) + MODIFICATION LIKE LINK GIVEN ABOVE

ecg_f = (reshape(ecg_f,1,length(ecg_f)));
[~,locs_R,delay]=pan_tompkin2(ecg_f,fs); %more accurate.
search_R = round(0.09*fs); % 90ms windowing time;
locs_Rn = locs_R - round(delay);

check_loc2 = find(locs_Rn>=length(ecg_f) | locs_Rn - search_R <=0 | locs_Rn + search_R >=length(ecg_f));
locs_Rn((check_loc2)) = [];
locs_R1 = locs_Rn + search_R;
locs_R2 = locs_Rn - search_R;


locs_Rf = zeros(1,length(locs_Rn));
for k = 1:length(locs_Rn)
    %%%%%%%% OLD: bug if there are two time points with same value
    % locs_Rf(k) = find(ecg_f == max(ecg_f(locs_R2(k):locs_R1(k))));
    %%%%%%%%%% PIERRE
    % Find the candidates for the peak
    temp_PC = find(ecg_f == max(ecg_f(locs_R2(k):locs_R1(k))));
    % If several values, find the one closer to the P peak
    if length(temp_PC) > 1
        % Compute dist to the P peak
        dist_with_R_peak = [];
        for iter_PC = 1:length(temp_PC)
            value_iter_PC = temp_PC(iter_PC);
            dist_with_R_peak = [dist_with_R_peak, abs(locs_R2(k) - value_iter_PC)];
        end
        % Find the min dist locs_Rf
        [minValue_PC, minIndex_PC] = min(dist_with_R_peak);
        % Add it to locs_Pf
        locs_Rf(k) = temp_PC(minIndex_PC);
    else
        locs_Rf(k) =  temp_PC;
    end
end
locs_Rf = unique(locs_Rf);
%% SEARCHING FOR Q PEAKS BEHIND R PEAKS. (MIN).
locs_Rf1 = locs_Rf;
search_offsetQ =  round(0.1*fs);%round(0.005*fs); %6
check_loc0 = find(locs_Rf1>=length(ecg_f) | locs_Rf1 - search_offsetQ <=0);
locs_Rf1((check_loc0)) = [];

locs_Qfround1 = locs_Rf1 - search_offsetQ; % searching offset.6
locs_Q = zeros(1,length(locs_Rf1));


for k = 1:length(locs_Rf1)
    [~,locs_Q(k)] =  min(ecg_f(locs_Qfround1(k):locs_Rf1(k)));
    locs_Q(k) = locs_Q(k) + locs_Qfround1(k);
end 
locs_Q = unique(locs_Q);
%% SEARCHING AHEAD OF R PEAKS FOR S PEAKS.
locs_Rf2 = locs_Rf;
search_offsetS = round(0.1*fs);
check_loc1 =  find(locs_Rf2>=length(ecg_f) | locs_Rf2 + search_offsetS >=length(ecg_f));
locs_Rf2((check_loc1)) = [];
locs_Sfround1 = locs_Rf2 + search_offsetS;

locs_Sf = zeros(1,length(locs_Rf2));
for k = 1:length(locs_Rf2)
    %%%%%%%% OLD: bug if there are two time points with same value
    % locs_Sf(k) =  find(ecg_f == min(ecg_f(locs_Rf2(k):locs_Sfround1(k))));
    %%%%%%%%%% PIERRE
    % Find the candidates for the peak
    temp_PC = find(ecg_f == min(ecg_f(locs_Rf2(k):locs_Sfround1(k))));
    % If several values, find the one closer to the P peak
    if length(temp_PC) > 1
        % Compute dist to the P peak
        dist_with_R_peak = [];
        for iter_PC = 1:length(temp_PC)
            value_iter_PC = temp_PC(iter_PC);
            dist_with_R_peak = [dist_with_R_peak, abs(locs_Rf2(k) - value_iter_PC)];
        end
        % Find the min dist locs_Rf2
        [minValue_PC, minIndex_PC] = min(dist_with_R_peak);
        % Add it to locs_Rf2
        locs_Sf(k) = temp_PC(minIndex_PC);
    else
        locs_Sf(k) =  temp_PC;
    end
end
locs_Sf = unique(locs_Sf);
%% SEARCHING FOR T PEAKS AFTER S PEAKS
locs_Tf1 = locs_Sf;
if length(locs_Rf)>1
    search_offsetT = round(mean(diff(locs_Rf))/2);
else
    search_offsetT = round(0.1*fs);
end 
check_loc2 = find(locs_Tf1>=length(ecg_f) | locs_Tf1 + search_offsetT >=length(ecg_f));
flag_T = 0;
locs_Tf = zeros(1,length(locs_Tf1));
if ~isempty(check_loc2)
    for k = 1:length(check_loc2)
        locs_Tf(check_loc2(k)) = find(ecg_f == max(ecg_f((locs_Tf1(check_loc2(k))) : length(ecg_f))),1);
    end
    flag_T = 1;
end 
if flag_T == 0
    locs_Tfround2 = locs_Tf1 + search_offsetT;
    locsfinale2 = find(locs_Tfround2>length(ecg_f));
    locs_Tfround2((locsfinale2)) = length(ecg_f);
    
    for k = 1:length(locs_Tf1)
        locs_Tf(k) =  find(ecg_f == max(ecg_f(locs_Tf1(k):locs_Tfround2(k))));
    end
elseif flag_T == 1
    locs_Tfround2 = locs_Tf1 + search_offsetT;
    locsfinale2 = find(locs_Tfround2>length(ecg_f));
    locs_Tfround2((locsfinale2)) = length(ecg_f);
    
    for k = 1:length(locs_Tf1)-length(check_loc2)
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        %%%%%%%%%% OLD (bug if two timepoint have exactly the same value
        % locs_Tf(k) =  find(ecg_f == max(ecg_f(locs_Tf1(k):locs_Tfround2(k))));
        %%%%%%%%%% PIERRE
        % Find the candidates for the peak
        temp_PC =  find(ecg_f == max(ecg_f(locs_Tf1(k):locs_Tfround2(k))));
        % If several values, find the one closer to the P peak
        if length(temp_PC) > 1 
            % Compute dist to the P peak
            dist_with_P_peak = [];
            for iter_PC = 1:length(temp_PC)
                value_iter_PC = temp_PC(iter_PC);
                dist_with_P_peak = [dist_with_P_peak, abs(locs_Tf1(k) - value_iter_PC)];
            end
            % Find the min dist locs_Pf
            [minValue_PC, minIndex_PC] = min(dist_with_P_peak);
            % Add it to locs_Pf
            locs_Tf(k) = temp_PC(minIndex_PC);
        else
            locs_Tf(k) =  temp_PC;
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    end
end 
locs_Tf = unique(locs_Tf);
%% SEARCHING FOR P PEAKS BEFORE Q PEAKS
locs_Pf1 = locs_Q;
locs_Pf = zeros(1,length(locs_Pf1));
if length(locs_Rf)>1
    search_offsetP = round(mean(diff(locs_Rf))/3);
else
    search_offsetP = round(0.1*fs);
end 
check_loc3 = find(locs_Pf1 - search_offsetP <=0);
flag_P = 0;
if ~isempty(check_loc3)
    for k = 1:length(check_loc3)
        locs_Pf(check_loc3(k)) = find(ecg_f == max(ecg_f((1:locs_Pf1(check_loc3(k))))));
    end 
        flag_P = 1;
end 
if flag_P == 0
    locs_Pfround1 = locs_Pf1 - search_offsetP;
    locsfinale = find(locs_Pfround1<=0);
    locs_Pfround1((locsfinale)) = 1;
    
    for k = 1:length(locs_Pf1)
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        %%%%%%%%%% OLD (bug if two timepoint have exactly the same value
        % locs_Pf(k) =  find(ecg_f == max(ecg_f(locs_Pfround1(k):locs_Pf1(k))));
        %%%%%%%%%% PIERRE
        % Find the candidates for the peak
        temp_PC = find(ecg_f == max(ecg_f(locs_Pfround1(k):locs_Pf1(k))));
        % If several values, find the one closer to the P peak
        if length(temp_PC) > 1 
            % Compute dist to the P peak
            dist_with_P_peak = [];
            for iter_PC = 1:length(temp_PC)
                value_iter_PC = temp_PC(iter_PC);
                dist_with_P_peak = [dist_with_P_peak, abs(locs_Pfround1(k) - value_iter_PC)];
            end
            % Find the min dist locs_Pf
            [minValue_PC, minIndex_PC] = min(dist_with_P_peak);
            % Add it to locs_Pf
            locs_Pf(k) = temp_PC(minIndex_PC);
        else
            locs_Pf(k) =  temp_PC;
        end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


    end
elseif flag_P == 1
    locs_Pfround1 = (locs_Pf1 - search_offsetP); 
    locsfinale = find(locs_Pfround1<=0);
    locs_Pfround1((locsfinale)) = 1;
    
    for k = 2:length(locs_Pf1)
        %%%%%%%% OLD: bug if there are two time points with same value
        % locs_Pf(k) =  find(ecg_f == max(ecg_f(locs_Pfround1(k):locs_Pf1(k))));
        %%%%%%%%%% PIERRE
        % Find the candidates for the peak
        temp_PC = find(ecg_f == max(ecg_f(locs_Pfround1(k):locs_Pf1(k))));
        % If several values, find the one closer to the P peak
        if length(temp_PC) > 1
            % Compute dist to the P peak
            dist_with_peak = [];
            for iter_PC = 1:length(temp_PC)
                value_iter_PC = temp_PC(iter_PC);
                dist_with_peak = [dist_with_peak, abs(locs_Pfround1(k) - value_iter_PC)];
            end
            % Find the min dist locs_Pf
            [minValue_PC, minIndex_PC] = min(dist_with_peak);
            % Add it to locs_Pf
            locs_Pf(k) = temp_PC(minIndex_PC);
        else
            locs_Pf(k) =  temp_PC;
        end
    end
end  
%% ADAPTIVE THRESHOLDING TO CONFIRM THAT THE PEAKS ARE INDEED THE P PEAKS (AVOIDING VERY SMALL P PEAKS IN AF SIGNALS)- flutter& af
T2 = zeros(1,min([length(locs_Pf),length(locs_Rf)]));
T1 = zeros(1,min(length(locs_Pf),length(locs_Rf)));
flag_Pfinal = zeros(1,min(length(locs_Pf),length(locs_Rf)));

for k1 = 1:(length(locs_Rf))
    T1(k1) = ecg_f(locs_Rf(k1))/10; % 2.5/0.25 mV as per standard values. (min)
    T2(k1) = ecg_f(locs_Rf(k1))/4; % trial and error(max).
end 

for k = 1:min([length(locs_Pf),length(locs_Rf),length(T1),length(T2)]) % for varying baseline we need more threshold
    if  (ecg_f(locs_Pf(k)) > T2(k)) 
        T2(k) = T2(k)*3.3;
    elseif ecg_f(locs_Pf(k)) <T1(k)
        T1(k) = T1(k)/2; % keep same.(2)
    end
end    
for k = 1:min([length(locs_Pf),length(locs_Rf),length(T1),length(T2)])
      if (ecg_f(locs_Pf(k)) < T1(k)) || (ecg_f(locs_Pf(k)) > T2(k))
          flag_Pfinal(k) = 1;
      end
end  
if plot_data == 1
    T1=interp1(1:length(T1),T1,linspace(1,length(T1),length(ecg_f)));
    T2=interp1(1:length(T2),T2,linspace(1,length(T2),length(ecg_f)));
end 

locs_to_eliminate = find(flag_Pfinal == 1);
locs_Pf((locs_to_eliminate)) = [];
locs_Pf = unique(locs_Pf);
%% CHECK IF R AND T PEAKS ARE OVERLAPPING IF SO ELIMINATE THE R PEAKS
[~,com,~] = intersect(locs_Rf,locs_Tf);
locs_Rf((com)) = [];
%% CHECKING FLAG ELEMENT TO MAKE SURE TO GET PLOTTING DATA FOR THRESHOLD.
if plot_data == 1
    varargout{1} = T1;
    varargout{2} = T2;
else
    varargout{1} = [];
    varargout{2} = [];
end
end 

