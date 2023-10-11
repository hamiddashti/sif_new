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
for k = 1:1:4
    for i = 1:length(leaf_ids)
        data = prepare_data(observed,"../data/","vilfan",leaf_ids(i));
        sample = data.observed_lrc;
        n = length(sample.Ca);

        obs.Ci = sample.Ci;
        obs.Ca = sample.Ca;
        obs.An = sample.An;
        obs.Fs = sample.Fs;
        obs.Fm_ = sample.Fm_;
        obs.Fm = sample.Fm;
        obs.Q = sample.PARi;
        obs.T = sample.T;

        obs.PhiP = 1- obs.Fs./obs.Fm_;  %PAM1
        obs.PhiN = obs.Fs.*(1./obs.Fm_-1./obs.Fm); %PAM2
        obs.PhiDF = obs.Fs./obs.Fm; %PAM3
        %% Set up inversion
        data.Q = obs.Q;
        data.T = obs.T;
        data.Ci = obs.Ci;
        data.O = repmat(209,n,1);
        data.E = 1e-03;
        v_inv = configure_fun_static(data);
        v_inv.Abs = 0.85;           % Total leaf absorptance to PAR, mol mol-1
        v_inv.beta = 0.52;          % PSII fraction of total absorptance, mol mol-1
        v_inv.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
        v_inv.CB6F = 175./v_inv.kq.*1e-06;      % Cyt b6f density, mol sites m-2
        v_inv.RUB = 50./v_inv.kc.*1e-06;      % Rubisco density, mol sites m-2
        v_inv.ss = solve_c3();

        % Set dynamic cross-sections and assign optimization function to 'v'
        v_inv.alpha_opt = 'static' ;
        % steps required for dynamic solver (doesnt need to be changed)
        if strcmp(v_inv.alpha_opt,'dynamic')==1
            v_inv.dynamic_solver = dynamic_solver();
            v_inv.c_steps = linspace(0,400e-06,10000);
            v_inv.gtc=repmat(0.05,n,1);
        end
        % Initial values
        e1e2_0=1;
        Vqmax_0 = 175;
        Vcmax_0 = 100;
        rm_0 = 12;
        vars_0 = [e1e2_0,Vqmax_0,Vcmax_0,rm_0]';

        % Range of values
        xmin = [0 20 10 1]';
        xmax = [5 500 200 12]';
        % What goes (==1) into objective function evaluation
        % CMAES algorithm settings
        opts.Restarts=3;
        opts.Noise.on=1;
        opts.LBounds = xmin;
        opts.UBounds= xmax;

        if k == 1
            weights.wAn = 1;
            weights.wPhiP = 0;
            weights.wPhiN =0;
            weights.wPhiDF = 0;
            [est_lrc_an(:,i),fval_cmaes_an(:,i)] = cmaes('chi2_static_real_lrc',vars_0,[],...
                opts,obs,weights,v_inv)
            
        elseif k ==2
            weights.wAn = 1;
            weights.wPhiP = 1;
            weights.wPhiN =0;
            weights.wPhiDF = 0;
            [est_lrc_an_phip(:,i),fval_cmaes_an_phip(:,i)] = cmaes('chi2_static_real_lrc',vars_0,[],...
                opts,obs,weights,v_inv)
            
        elseif k ==3
            weights.wAn = 1;
            weights.wPhiP = 1;
            weights.wPhiN =1;
            weights.wPhiDF = 0;
            [est_lrc_an_phip_phin(:,i),fval_cmaes_an_phip_phin(:,i)] = cmaes('chi2_static_real_lrc',vars_0,[],...
                opts,obs,weights,v_inv)
            
        elseif k==4
            weights.wAn = 1;
            weights.wPhiP = 1;
            weights.wPhiN =1;
            weights.wPhiDF = 1;
            [est_lrc_an_phip_phidf_phidf(:,i),fval_cmaes_an_phip_phin_phidf(:,i)] = cmaes('chi2_static_real_lrc',vars_0,[],...
                opts,obs,weights,v_inv)
            
        end

    end
end
save('est_lrc_an.mat','est_lrc_an')
save('est_lrc_an_phip.mat','est_lrc_an_phip')
save('est_lrc_an_phip_phin.mat','est_lrc_an_phip_phin')
save('est_lrc_an_phip_phidf_phidf.mat','est_lrc_an_phip_phidf_phidf')
