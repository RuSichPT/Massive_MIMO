clc; clear;%close all
Eb_N0 = 0:40;
%%
figure(1)
load('DataBase/4x4x4_pre =0_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/ScatteringFlat/8x4x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/ScatteringFlat/16x4x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/ScatteringFlat/32x4x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
load('DataBase/ScatteringFlat/64x4x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-',1.5,0,[0.45 0.45 0.45]);
legend('4x4x4','8x4x4','16x4x4','32x4x4','64x4x4');
ylim([10^-5 10^0]);
%%
figure(2)
load('DataBase/4x4x4_pre =0_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/ScatteringFlat/4x8x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/ScatteringFlat/4x16x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/ScatteringFlat/4x32x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
load('DataBase/ScatteringFlat/4x64x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-',1.5,0,[0.45 0.45 0.45]);
legend('4x4x4','4x8x4','4x16x4','4x32x4','4x64x4');
ylim([10^-5 10^0]);
%%
figure(3)
load('DataBase/4x4x4_pre =0_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/8x8x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/16x16x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/32x32x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
load('DataBase/64x64x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-',1.5,0,[0.45 0.45 0.45]);
legend('4x4x4','8x8x4','16x16x4','32x32x4','64x64x4');
ylim([10^-5 10^0]);
%%
figure(4)
load('DataBase/4x4x4_pre =0_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/16x8x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/32x16x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/64x32x4_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
legend('4x4x4','16x8x4','32x16x4','64x32x4');
ylim([10^-5 10^0]);
%%
figure(5)
load('DataBase/8x8x8_pre =0_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'k',1.5,0);
load('DataBase/32x16x8_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'--k',1.5,0);
load('DataBase/64x16x8_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,'-.k',1.5,0);
load('DataBase/64x32x8_pre =1_ster =0_ScatteringFlat_Wmm=0_Wm=0_Exp=1.mat');
mean_bear = mean(ber_mean,1);
plot_ber(mean_bear,SNR,prm.bps*prm.numSTS,':k',1.5,0);
legend('8x8x8','32x16x8','64x16x8','64x32x8');
ylim([10^-5 10^0]);