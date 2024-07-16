function z = baseline_remove(sig)
%% what does this function do?
% This is a DC notch filter to remove baseline noise from ECG signals for
% further analysis
%% logic
b = [1 -1];
a = [1 -0.9];
z = filter(b,a,sig);
end 
