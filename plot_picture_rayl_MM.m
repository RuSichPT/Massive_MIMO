clc; clear;%close all
Eb_N0 = 0:40;
%% ЧБ
%% beamfrming и передача копий
figure(1)
load('DataBase/RAYL/4x4x4_pre =0_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL/8x4x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/RAYL/16x4x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/RAYL/32x4x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
legend('4x4x4','8x4x4','16x4x4','32x4x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% разнесенный прием
figure(2)
load('DataBase/RAYL/4x4x4_pre =0_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL/4x8x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/RAYL/4x16x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/RAYL/4x32x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
legend('4x4x4','4x8x4','4x16x4','4x32x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% симметричные антенны
figure(3)
load('DataBase/RAYL/4x4x4_pre =0_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL/8x8x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/RAYL/16x16x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/RAYL/32x32x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
legend('4x4x4','8x8x4','16x16x4','32x32x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% ЦВЕТ
%% beamfrming и передача копий
figure(4)
load('DataBase/RAYL/4x4x4_pre =0_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL/8x4x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'r',1.5,0);
load('DataBase/RAYL/16x4x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'b',1.5,0);
load('DataBase/RAYL/32x4x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'g',1.5,0);
legend('4x4x4','8x4x4','16x4x4','32x4x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% разнесенный прием
figure(5)
load('DataBase/RAYL/4x4x4_pre =0_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL/4x8x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'r',1.5,0);
load('DataBase/RAYL/4x16x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'b',1.5,0);
load('DataBase/RAYL/4x32x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'g',1.5,0);
legend('4x4x4','4x8x4','4x16x4','4x32x4');
ylim([10^-5 10^0]);
xlim([0 40]);
%% симметричные антенны
figure(6)
load('DataBase/RAYL/4x4x4_pre =0_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/RAYL/8x8x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'r',1.5,0);
load('DataBase/RAYL/16x16x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'b',1.5,0);
load('DataBase/RAYL/32x32x4_pre =1_ster =0_RAYL_Wmm=0_Wm=0_Exp=50.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'g',1.5,0);
legend('4x4x4','8x8x4','16x16x4','32x32x4');
ylim([10^-5 10^0]);
xlim([0 40]);