% Clean up working environment
clear all; % variables
close all; % figures
clc; % command window

homedir = '/home/hamid/SIF/codes/sif_new/';
outdir = strcat(homedir,"outputs/");
addpath(genpath(strcat(homedir,'src/')));


%% Select inputs

% Specify environmental conditions
n = 10;                                   % Steps in vector
data.Q = repmat(1200,n,1);   % PAR, umol PPFD m-2 s-1
data.T = repmat(25,n,1);                  % Leaf temperature, C
data.Ci = transpose(linspace(50,1000,n));                 % Atmospheric CO2, ubar
data.O = repmat(209,n,1);                 % Atmospheric O2, mbar
data.E = 1e-03;

%% Configure model simulations

% Load default set of parameters
v_obs = configure_fun_static(data);

% Adjust parameters of interest from default values
v_obs.Abs = 0.85;           % Total leaf absorptance to PAR, mol mol-1
v_obs.beta = 0.52;          % PSII fraction of total absorptance, mol mol-1
v_obs.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
v_obs.CB6F = 175./v_obs.kq.*1e-06;      % Cyt b6f density, mol sites m-2
v_obs.RUB = 50./v_obs.kc.*1e-06;      % Rubisco density, mol sites m-2
v_obs.ss = solve_c3();

v_obs.rm = 0.01;
v_obs.e1e2 = 2;
v_obs.eps1 = v_obs.e1e2./(1+v_obs.e1e2); % PS I transfer function, mol mol-1
v_obs.eps2 = 1./(1+v_obs.e1e2);
v_obs.Vqmax = 175;
v_obs.Vcmax = 100;

%% Run simulation and visualize results

% Set dynamic cross-sections and assign optimization function to 'v' 
v_obs.alpha_opt = 'static' ;
% steps required for dynamic solver (doesnt need to be changed) 
if strcmp(v_obs.alpha_opt,'dynamic')==1
v_obs.dynamic_solver = dynamic_solver();
v_obs.c_steps = linspace(0,400e-06,10000);
v_obs.gtc=repmat(0.05,n,1);
end


%%
% Run simulation
observed_static_co2 = model_fun_static(v_obs);
%%
figure(2)
subplot(5,1,1)
plot(data.Ci,observed_static_co2.An_a*1e6,'*-b')
ylabel('An_a',FontSize=13,FontWeight='bold')
xlabel('PAR',FontSize=13,FontWeight='bold')
subplot(5,1,2)
plot(data.Ci,observed_static_co2.PAM1_a,'*-b')
ylabel('PhiP',FontSize=13,FontWeight='bold')
xlabel('PAR',FontSize=13,FontWeight='bold')
subplot(5,1,3)
plot(data.Ci,observed_static_co2.PAM2_a,'*-b')
ylabel('PhiN',FontSize=13,FontWeight='bold')
xlabel('PAR',FontSize=13,FontWeight='bold')
subplot(5,1,4)
plot(data.Ci,observed_static_co2.PAM3_a,'*-b')
ylabel('PhiD+PhiF',FontSize=13,FontWeight='bold')
xlabel('PAR',FontSize=13,FontWeight='bold')

subplot(5,1,5)
plot(data.Ci,observed_static_co2.a2./v_obs.Abs.*100,'*-b')
ylabel('a2',FontSize=13,FontWeight='bold')
xlabel('PAR',FontSize=13,FontWeight='bold')
hold on 
plot(data.Ci,observed_static_co2.a1./v_obs.Abs.*100,'*-r')
ylabel('a1',FontSize=13,FontWeight='bold')
xlabel('PAR',FontSize=13,FontWeight='bold')
ylim([0 100]);

% a2 = observed_syn.a2;
% a1 =  observed_syn.a1;
% C = observed_syn.C;
observed_static_co2.v_obs = v_obs;
save(strcat(outdir,"observed_static_co2.mat"),"observed_static_co2")
% save(strcat(homedir,"C.mat"),"C")