clc
clear all
close all

fs = 360;
try
    ecg_1=load('/Users/rohansanghavi/Desktop/ECG_PQRST_VERSION_3/8 Bigeminy/106m (8).mat');
catch 
    error('File is not found');
end 
plot_data = 1; % (REFER FUNCTION DOCUMENTATION)
ecg_1 = ecg_1.val;
ecg_f = preprocess_window_ecg(ecg_1,fs);

[locs_P,locs_Q,locs_R,locs_S,locs_T,T1,T2] = compute_fudicial_peaks_live17_c(ecg_f,fs,plot_data);

[mean_QR,mean_SS,mean_ST,mean_RS,std_QR,std_SS,std_ST,skew_QR,skew_SS,skew_ST,skew_RS,kurt_QR,kurt_SS,kurt_ST,kurt_RS, std_RS] = compute_mean_intervals2_c(locs_R,locs_S,locs_Q,locs_T,fs);

[RR_mean,QRS_mean,mean_QT,mean_PR,RR_std,RR_skew,RR_kurt,QRS_std,QRS_skew,QRS_kurt,std_PR,skew_PR,kurt_PR,std_QT,skew_QT,kurt_QT]= compute_mean_interval_c(locs_P,locs_Q,locs_R,locs_S,locs_T,fs);

bpm = round(60/RR_mean);


figure(1);
hold on
plot(ecg_f,'linewidth',1.5);
plot(locs_R,ecg_f(locs_R),'rv','MarkerFaceColor','r');
plot(locs_Q,ecg_f(locs_Q),'rs','MarkerFaceColor','g');
plot(locs_S,ecg_f(locs_S),'o','MarkerFaceColor','y');
plot(locs_P,ecg_f(locs_P),'v','MarkerFaceColor','m');
plot(locs_T,ecg_f(locs_T),'s','MarkerFaceColor','k');
plot(T2,'linewidth',2,'Linestyle','--','Color','red');
plot(T1,'linewidth',2,'Linestyle','--','Color','green');
legend('ECG Signal','R peaks','Q Peaks','S Peaks','P peaks','T peaks','T-P1','T-P2');
hold off
title('QRS complex being marked on the fitered ECG signal along with the P and T peaks');
xlabel('Samples');
ylabel('Normalised Amplitude');
grid on




