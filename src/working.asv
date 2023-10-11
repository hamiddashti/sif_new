clear all
clc
close all
homedir = '/home/hamid/SIF/codes/sif_new/';
outdir = strcat(homedir,"outputs/");
addpath(genpath(strcat(homedir)));

% Read Vilfan's co2 and lrc data
% observed.vilfan_data = prepare_vilfan_data("../data/");
% save("observed_vilfan.mat","observed")
load("observed_vilfan.mat")
leaf_ids = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
%% Load the LRC observation

data = prepare_data(observed,"../data/","vilfan",leaf_ids(7));
sample = data.observed_co2;
n = length(sample.Ca);

obs.Ci = sample.Ci;
obs.Ca = sample.Ca;
obs.An = sample.An;
obs.Fs = sample.Fs;
obs.Fm_ = sample.Fm_;
obs.Fm = sample.Fm;
obs.Q = sample.PARi;
obs.T = sample.T;
obs.ETR = sample.ETR

obs.PhiP = 1- obs.Fs./obs.Fm_;  %PAM1
obs.PhiN = obs.Fs.*(1./obs.Fm_-1./obs.Fm); %PAM2
obs.PhiDF = obs.Fs./obs.Fm; %PAM3
subplot(2,1,1)
plot(obs.Ci,obs.An)
subplot(2,1,2)
plot(obs.Ci,obs.ETR)