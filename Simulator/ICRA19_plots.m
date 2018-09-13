
clear all; close all; clc;

start_epoch= 10;
end_epoch= 220;

alert_limit= 1;

% load DATA
load('DATA/M2_s30_x15_y4.mat');
P_HMI_M2= DATA.P_HMI(start_epoch:end_epoch);

load('DATA/M3_s30_x15_y4.mat');
P_HMI_M3= DATA.P_HMI(start_epoch:end_epoch);

load('DATA/M4_s30_x15_y4.mat');
P_HMI_M4= DATA.P_HMI(start_epoch:end_epoch);

load('DATA/M5_s30_x15_y4.mat');
P_HMI_M5= DATA.P_HMI(start_epoch:end_epoch);

load('DATA/M6_s30_x15_y4.mat');
P_HMI_M6= DATA.P_HMI(start_epoch:end_epoch);



% find P_HMI_H0
for i= start_epoch:end_epoch
    P_HMI_H0(i)= 2* normcdf(-alert_limit,0,DATA.stdEps(i)/3);
end
P_HMI_H0= P_HMI_H0(start_epoch:end_epoch);



% Plots
figure; hold on; grid on;
set(gca,'Yscale','log');
plot(start_epoch:end_epoch, P_HMI_M2, 'linewidth', 2)
% plot(start_epoch:end_epoch, P_HMI_M3, 'linewidth', 2)
plot(start_epoch:end_epoch, P_HMI_M4, 'linewidth', 2)
% plot(start_epoch:end_epoch, P_HMI_M5, 'linewidth', 2)
plot(start_epoch:end_epoch, P_HMI_M6, 'linewidth', 2)
% plot(start_epoch:end_epoch, P_HMI_H0, 'linewidth', 2)
% legend('M = 2', 'M = 3', 'M = 4', 'M = 5', 'M = 6') 
legend('M = 2','M = 4', 'M = 6') 

