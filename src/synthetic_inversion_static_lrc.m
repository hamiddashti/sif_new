%% Clean up working environment
clear all; % variables
close all; % figures
clc; % command window

homedir = '/home/hamid/SIF/codes/sif_new/';
outdir = strcat(homedir,"outputs/");
addpath(genpath(strcat(homedir)))
load(strcat(outdir,'observed_static.mat'))
% load(strcat(homedir,'a2.mat'))
% load(strcat(homedir,'a1.mat'))
% load(strcat(homedir,'C.mat'))

%% Set up environment
n = 15;                                   % Steps in vector
data.Q = transpose(linspace(10,2400,n));   % PAR, umol PPFD m-2 s-1
data.T = repmat(25,n,1);                  % Leaf temperature, C
% data.Ca = repmat(400,n,1);                 % Atmospheric CO2, ubar
data.Ci = repmat(270,n,1);
data.O = repmat(209,n,1);                 % Atmospheric O2, mbar
data.E = 1e-03;

%% Set other paramteres 

% Load default set of parameters
v_inv = configure_fun_static(data);

% Adjust parameters of interest from default values
v_inv.Abs =observed_static.v_obs.Abs;           % Total leaf absorptance to PAR, mol mol-1
v_inv.beta = observed_static.v_obs.beta;          % PSII fraction of total absorptance, mol mol-1
v_inv.Ku2 = observed_static.v_obs.Ku2;           % Rate constant for exciton sharing at PSII, s-1
v_inv.CB6F = observed_static.v_obs.CB6F;      % Cyt b6f density, mol sites m-2
v_inv.RUB = observed_static.v_obs.RUB;      % Rubisco density, mol sites m-2
v_inv.ss = solve_c3();
% v_inv.rm = 8;
% v_inv.e1e2 = 2;
% v_inv.eps1 = v_inv.e1e2./(1+v_inv.e1e2); % PS I transfer function, mol mol-1
% v_inv.eps2 = 1./(1+v_inv.e1e2);
% v_inv.Vqmax = 175;
% v_inv.Vcmax = 100;

%%
v_inv.alpha_opt = 'static';
% steps required for dynamic solver (doesnt need to be changed) 
if strcmp(v_inv.alpha_opt,'dynamic')==1
v_inv.dynamic_solver = dynamic_solver();
v_inv.c_steps = linspace(0,400e-06,10000);
v_inv.gtc=repmat(0.05,n,1);
end
% if strcmp(v_inv.alpha_opt,'free')==1
%     v_inv.a2 = a2;
%     v_inv.a1 = a1;
% end
% 
%% Invert the model
% X_labels = {'e1e2','Vqmax','Vcmax','rm'};
% M=length(X_labels);
e1e2_0=2;
Vqmax_0 = 175;
Vcmax_0 = 100;
rm_0 = 8;

vars_0 = [e1e2_0,Vqmax_0,Vcmax_0,rm_0]';
xmin = [0 20 10 4]';
xmax = [5 500 200 12]';
weights.wAn = 1;
weights.wPhiP = 0;
weights.wPhiN =0;
weights.wPhiDF = 0;
weights.wETR = 0;

true_pars = [observed_static.v_obs.e1e2,observed_static.v_obs.Vqmax,...
    observed_static.v_obs.Vcmax, observed_static.v_obs.rm]'

%% GA algorithm
fh = @(vars)chi2_static(vars,observed_static,weights,v_inv);
% options = optimoptions('ga','Display','iter','UseParallel',true,...
%     'PopulationSize',1000,...
%     'InitialPopulationRange', [xmin(1) xmax(1); xmin(2) xmax(2);...
%     xmin(3) xmax(3);xmin(4) xmax(4)]');
% options = optimoptions('ga','Display','iter','UseParallel',true,...
%     'PopulationSize',500,'SelectionFcn','selectiontournament',...
%     'FitnessScalingFcn',@fitscalingprop);
% 
options = optimoptions('ga','Display','iter','UseParallel',true,...
    'PopulationSize',5000);

thestate = rng;
[est_ga,fval_ga,exitflag,output] = ga(fh,length(vars_0),[],[],[],[],xmin,...
    xmax,[],options)

[est_ga', true_pars]
% plot(est_ga(5:end)',true_pars(5:end),'*')
% 100*abs((est_ga'-true_pars))./true_pars


%% CMAES algorithm
opts.Restarts=3;
opts.Noise.on=1;
opts.LBounds = xmin;
opts.UBounds= xmax;

[est_cmaes,fval_cmaes] = cmaes('chi2_static',vars_0,[],...
    opts,observed_static,weights,v_inv)

[est_cmaes, true_pars]
% 100*abs((est_cmaes-true_pars))./true_pars
