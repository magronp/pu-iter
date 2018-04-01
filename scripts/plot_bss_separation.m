% Plot results of source separation benchmark
clear all; close all; clc;
test_or_dev = 'Test';
set_settings_puiter;

% Load data
load(strcat(metrics_path,'score_bss_oracle.mat')); SDRo=squeeze(mean(SDR,1))'; SIRo=squeeze(mean(SIR,1))'; SARo=squeeze(mean(SAR,1))';
load(strcat(metrics_path,'score_bss_informed.mat')); SDRi=squeeze(mean(SDR,1))'; SIRi=squeeze(mean(SIR,1))'; SARi=squeeze(mean(SAR,1))';
load(strcat(metrics_path,'score_bss_blind.mat')); SDRb=squeeze(mean(SDR,1))'; SIRb=squeeze(mean(SIR,1))'; SARb=squeeze(mean(SAR,1))';

% Iterative method
figure;
h11=subplot(3,3,1); boxplot(SDRo); title('SDR (dB)'); set(gca,'FontSize',14,'XtickLabel',[]);
h12=subplot(3,3,2); boxplot(SIRo); title('SIR (dB)'); set(gca,'FontSize',14,'XtickLabel',[]);
h13=subplot(3,3,3); boxplot(SARo); title('SAR (dB)'); set(gca,'FontSize',14,'XtickLabel',[]);

h21=subplot(3,3,4); boxplot(SDRi); set(gca,'FontSize',14,'XtickLabel',[]);
h22=subplot(3,3,5); boxplot(SIRi); set(gca,'FontSize',14,'XtickLabel',[]);
h23=subplot(3,3,6); boxplot(SARi); set(gca,'FontSize',14,'XtickLabel',[]);

h31=subplot(3,3,7); boxplot(SDRb); set(gca,'FontSize',14,'XtickLabel',algos,'XtickLabelRotation',90);
h32=subplot(3,3,8); boxplot(SIRb); set(gca,'FontSize',14,'XtickLabel',algos,'XtickLabelRotation',90);
h33=subplot(3,3,9); boxplot(SARb); set(gca,'FontSize',14,'XtickLabel',algos,'XtickLabelRotation',90);

boxheight = 0.22;

h11.Position(4) = boxheight; h12.Position(4) = boxheight; h13.Position(4) = boxheight;
h11.Position(2) = 0.7; h12.Position(2) = 0.7; h13.Position(2) = 0.7;

h21.Position(4) = boxheight; h22.Position(4) = boxheight; h23.Position(4) = boxheight;
h21.Position(2) = 0.45; h22.Position(2) = 0.45; h23.Position(2) = 0.45;

h31.Position(4) = boxheight; h32.Position(4) = boxheight; h33.Position(4) = boxheight;
h31.Position(2) = 0.2; h32.Position(2) = 0.2; h33.Position(2) = 0.2;
