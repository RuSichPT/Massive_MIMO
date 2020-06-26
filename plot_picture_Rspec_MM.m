clc; clear;%close all
Eb_N0 = 0:40;
%% beamfrming и передача копий
figure()
load('DataBase/RAYL_SPECIAL/4x4x4_pre =0_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL_SPECIAL/8x4x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/RAYL_SPECIAL/16x4x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/RAYL_SPECIAL/32x4x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
load('DataBase/RAYL_SPECIAL/64x4x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-',1.5,0,[0.45 0.45 0.45]);
legend('4x4x4','8x4x4','16x4x4','32x4x4','64x4x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% разнесенный прием
figure()
load('DataBase/RAYL_SPECIAL/4x4x4_pre =0_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL_SPECIAL/4x8x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/RAYL_SPECIAL/4x16x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/RAYL_SPECIAL/4x32x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
load('DataBase/RAYL_SPECIAL/4x64x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-',1.5,0,[0.45 0.45 0.45]);
legend('4x4x4','4x8x4','4x16x4','4x32x4','4x64x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% симметричные антенны
figure()
load('DataBase/RAYL_SPECIAL/4x4x4_pre =0_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL_SPECIAL/8x8x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/RAYL_SPECIAL/16x16x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/RAYL_SPECIAL/32x32x4_pre =1_ster =0_RAYL_SPECIAL_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
legend('4x4x4','8x8x4','16x16x4','32x32x4');
ylim([10^-5 10^0]);
xlim([0 40]);