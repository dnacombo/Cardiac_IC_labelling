function [RR_mean,QRS_mean,mean_QT,mean_PR,RR_std,RR_skew,RR_kurt,QRS_std,QRS_skew,QRS_kurt,std_PR,skew_PR,kurt_PR,std_QT,skew_QT,kurt_QT]= compute_mean_interval_c(locs_P,locs_Q,locs_R2,locs_S,locs_T,fs)
%% WHAT DOES THIS FUNCTION DO?
% THIS FUNCTION COMPUTES THE MEAN RR, MEAN QRS, MEAN QT AND MEAN PR
% INTERVAL BASED ON THE LOCATIONS OF THE PEAKS.
% OTHER STATISTICAL DATA ALSO COMPUTED.

% COMPILED AND MAINTAINED BY-
% ROHAN SANGHAVI
%% CODE LOGIC - RR INTERVAL
RR_period = diff((locs_R2))./fs;
RR_mean = mean(RR_period); % we consider the unwindowed mean for accuracy
RR_std =  std(RR_period);
RR_skew =  skewness(RR_period);
RR_kurt =  kurtosis(RR_period);

if isnan(RR_mean)
    RR_mean = 0;
end 

if isnan(RR_std)
    RR_std = 0;
end  

if isnan(RR_skew)
    RR_skew = 0;
end 

if isnan(RR_kurt)
    RR_kurt = 0;
end 
%% QRS COMPLEX PERIOD.
if  (length(unique(locs_Q))>2 && length(unique(locs_S))>2)
    if length(locs_Q)~=length(locs_S)
        locs_S = locs_S(~isnan(interp1(unique(locs_Q),unique(locs_Q),unique(locs_S),'nearest'))); % if lengths are unequal
    end
end

if  (length(unique(locs_Q))>2 && length(unique(locs_S))>2)
    if length(locs_Q)~=length(locs_S)
        locs_Q = locs_Q(~isnan(interp1(unique(locs_S),unique(locs_S),unique(locs_Q),'nearest'))); % interpolation to find closest values.
    end
end


% double check and then finally subtract.
if length(locs_Q) > length(locs_S)
    locs_Q1 = locs_Q(1:length(locs_S));
    locs_S1 = locs_S;
    
else 
    locs_S1 = locs_S(1:length(locs_Q));
    locs_Q1 = locs_Q;
   
end 
QRS_period = ((locs_S1 - locs_Q1)./fs);
QRS_mean = mean(QRS_period);
QRS_std = std(QRS_period);
QRS_kurt = kurtosis(QRS_period);
QRS_skew = skewness(QRS_period);


if isnan(QRS_mean)
    QRS_mean = 0;
end 

if isnan(QRS_std)
    QRS_std = 0;
end 

if isnan(QRS_skew)
    QRS_skew = 0;
end 

if isnan(QRS_kurt)
    QRS_kurt = 0;
end 
%% PR INTERVAL
if  (length(unique(locs_P))>2 && length(unique(locs_R2))>2)
    if length(locs_P)~=length(locs_R2)
        locs_R2 = locs_R2(~isnan(interp1(unique(locs_P),unique(locs_P),unique(locs_R2),'nearest'))); % if lengths are unequal
    end
end 
if  (length(unique(locs_P))>2 && length(unique(locs_R2))>2)
    if length(locs_P)~=length(locs_R2)
        locs_P = locs_P(~isnan(interp1(unique(locs_R2),unique(locs_R2),unique(locs_P),'nearest'))); % interpolation to find closest values.
    end
end 

if length(locs_R2) > length(locs_P)
    locs_R3 = locs_R2(1:length(locs_P));
    locs_P2 = locs_P;
    temp = 1;
else 
    locs_P2 = locs_P(1:length(locs_R2));
    locs_R3 = locs_R2;
    temp = 1;
end 

PR_interval = ((locs_R3 - locs_P2)./fs);
if temp == 0
    mean_PR = mean(PR_interval);
else
    mean_PR = abs(mean(PR_interval));
end 
std_PR =  std(PR_interval);
skew_PR =  skewness(PR_interval);
kurt_PR =  kurtosis(PR_interval);

if isnan(mean_PR)
    mean_PR = 0;
end 

if isnan(std_PR)
    std_PR = 0;
end

if isnan(kurt_PR)
    kurt_PR = 0;
end 

if isnan(skew_PR)
    skew_PR = 0;
end 

%% QT INTERVAL
if  (length(unique(locs_Q))>2 && length(unique(locs_T))>2)
    if length(locs_Q)~=length(locs_T)
        locs_T = locs_T(~isnan(interp1(unique(locs_Q),unique(locs_Q),unique(locs_T),'nearest'))); % if lengths are unequal
    end
end 

if  (length(unique(locs_Q))>2 && length(unique(locs_T))>2)
    if length(locs_Q)~=length(locs_T)
        locs_Q = locs_Q(~isnan(interp1(unique(locs_T),unique(locs_T),unique(locs_Q),'nearest'))); % interpolation to find closest values.
    end
end 

if length(locs_Q) > length(locs_T)
    locs_Q2 = locs_Q(1:length(locs_T));
    locs_T2 = locs_T;
    flag =1;
    
else 
    locs_T2 = locs_T(1:length(locs_Q));
    locs_Q2 = locs_Q;
    flag =1; % equal here
end 


QT_interval = ((locs_T2 - locs_Q2)./fs);
if flag ==0 
    mean_QT = mean(QT_interval);
elseif flag == 1
    mean_QT = abs(mean(QT_interval));
else
    mean_QT = abs(mean(QT_interval));
end 

std_QT = std(QT_interval);
skew_QT = skewness(QT_interval);
kurt_QT = kurtosis(QT_interval);

if isnan(mean_QT) 
    mean_QT = 0;
end 

if isnan(std_QT) 
    std_QT = 0;
end 

if isnan(skew_QT) 
    skew_QT = 0;
end 

if isnan(kurt_QT) 
    kurt_QT = 0;
end 
end 