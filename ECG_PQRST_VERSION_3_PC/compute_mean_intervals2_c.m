function [mean_QR,mean_SS,mean_ST,mean_RS,std_QR,std_SS,std_ST,skew_QR,skew_SS,skew_ST,skew_RS,kurt_QR,kurt_SS,kurt_ST,kurt_RS, std_RS] = compute_mean_intervals2_c(locs_R,locs_S,locs_Q,locs_T,fs)
%% WHAT DOES THIS FUNCTION DO?
% THIS FUNCTION COMPUTES THE MEAN QR,SS,ST AND RS INTERVALS IN AN ECG
% SIGNAL. 
% SS-> SIMPLY DIFFERENTIATE AND DIVIDE BY fs.

% QR -> WE USE THE arrayfun INBUILT FUNCTION IN MATLAB, TO FIND THE
% MINIMAL Q PEAK AFTER A SINGLE R PEAK. IN OTHER WORDS, IT ENSURES THAT
% THERE IS A MAPPING OF THE PEAKS.

% RS -> SAME LIKE QR INTERVAL

% ST-> WE JUST SIMPLY SUBTRACT AFTER TRUNCATING IN CASE THERE IS AN ERROR.
% ALSO STANDARD DEVIATION, SKEWNESS AND KURTOSIS ARE CALCULATED.
%% AUTHORS
% COMPILED AND MAINTAINED BY-
% ROHAN SANGHAVI.
%% SS INTERVAL CODE
SS_int = diff(locs_S)./fs;
mean_SS = mean(SS_int);
std_SS = std(SS_int);
skew_SS = skewness(SS_int);
kurt_SS = kurtosis(SS_int);
%% QR INTERVAL CODE
if  (length(unique(locs_Q))>2 && length(unique(locs_R))>2)
    if length(locs_Q)~=length(locs_R)
        locs_R = locs_R(~isnan(interp1(unique(locs_Q),unique(locs_Q),unique(locs_R),'nearest'))); % if lengths are unequal
    end
end 

if  (length(unique(locs_Q))>2 && length(unique(locs_R))>2)
    if length(locs_Q)~=length(locs_R)
        locs_Q = locs_Q(~isnan(interp1(unique(locs_R),unique(locs_R),unique(locs_Q),'nearest'))); % interpolation to find closest values.
    end
end 

if length(locs_Q) > length(locs_R)
    locs_Q = locs_Q(1:length(locs_R));
elseif length(locs_R)>length(locs_Q)
    locs_R = locs_R(1:length(locs_Q));
end 
QR_int = (locs_R - locs_Q)./fs;
mean_QR = mean(QR_int);
std_QR = std(QR_int);
skew_QR = skewness(QR_int);
kurt_QR = kurtosis(QR_int);
%% ST INTERVAL CODE
if  (length(unique(locs_S))>2 && length(unique(locs_T))>2)
    if length(locs_S)~=length(locs_T)
        locs_T = locs_T(~isnan(interp1(unique(locs_S),unique(locs_S),unique(locs_T),'nearest'))); % if lengths are unequal
    end
end

if  (length(unique(locs_S))>2 && length(unique(locs_T))>2)
    if length(locs_S)~=length(locs_T)
        locs_S = locs_S(~isnan(interp1(unique(locs_T),unique(locs_T),unique(locs_S),'nearest'))); % interpolation to find closest values.
    end
end

if length(locs_T)<length(locs_S)
    locs_Sf = locs_S(1:length(locs_T));
    locs_Tf = locs_T;
    
elseif length(locs_T)>length(locs_S)
    locs_Tf = locs_T(1:length(locs_S));
    locs_Sf = locs_S;
else
    locs_Tf = locs_T;
    locs_Sf = locs_S;
end 

ST_int = (locs_Tf - locs_Sf)./fs;
mean_ST = mean(ST_int);
std_ST = std(ST_int);
skew_ST  = skewness(ST_int);
kurt_ST = kurtosis(ST_int);
%% SR INTERVAL CODE
if  (length(unique(locs_S))>2 && length(unique(locs_R))>2)
    if length(locs_R)~=length(locs_S)
        locs_R = locs_R(~isnan(interp1(unique(locs_S),unique(locs_S),unique(locs_R),'nearest'))); % if lengths are unequal
    end
end 

if  (length(unique(locs_S))>2 && length(unique(locs_R))>2)
    if length(locs_R)~=length(locs_S)
        locs_S = locs_S(~isnan(interp1(unique(locs_R),unique(locs_R),unique(locs_S),'nearest'))); % interpolation to find closest values.
    end
end 

if length(locs_S) < length(locs_R)
    locs_R = locs_R(1:length(locs_S));
elseif length(locs_S) > length(locs_R)
    locs_S = locs_S(1:length(locs_R));
end 
RS_int = (locs_S - locs_R)./fs;
mean_RS = mean(RS_int);
std_RS = std(RS_int);
skew_RS = skewness(RS_int);
kurt_RS = kurtosis(RS_int);
%% ELIMINATING UNWANTED POSSIBILITIES.
if isnan(mean_QR) 
    mean_QR = 0;
end

if isnan(mean_SS) 
    mean_SS = 0;
end

if isnan(mean_ST) 
    mean_ST = 0;
end

if isnan(mean_RS) 
    mean_RS = 0;
end




if isnan(std_QR) 
    std_QR = 0;
end

if isnan(std_SS) 
    std_SS = 0;
end

if isnan(std_ST) 
    std_ST = 0;
end

if isnan(std_RS) 
    std_RS = 0;
end



if isnan(skew_QR) 
    skew_QR = 0;
end

if isnan(skew_SS) 
    skew_SS = 0;
end

if isnan(skew_ST) 
    skew_ST = 0;
end

if isnan(skew_RS) 
    skew_RS = 0;
end



if isnan(kurt_QR) 
    kurt_QR = 0;
end

if isnan(kurt_SS) 
    kurt_SS = 0;
end

if isnan(kurt_ST) 
    kurt_ST = 0;
end

if isnan(kurt_RS) 
    kurt_RS = 0;
end


end 







