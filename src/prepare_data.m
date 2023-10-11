function observed = prepare_data(observed_real,datadir,dataset,sample_number)

% load(strcat(homedir,"observed_real.mat"))
% load(strcat(homedir,"spectral.mat"))

if strcmp(dataset,'will')==1

    abs_lrc = readtable(strcat(datadir,'william_woodgate/LRC_leaf_abs.xlsx'),...
        'Sheet','LRC_abs');

    % Prepare CO2 experiment for the sample

    % Gas exchange data
    licor_co2 = observed_real.aus_data.licor_co2;
    spec_co2 = observed_real.aus_data.spec_co2;
    observed_co2.Ca = licor_co2.Ca(:,sample_number);
    observed_co2.Ci = licor_co2.Ci(:,sample_number);
    observed_co2.Q = licor_co2.Qin(:,sample_number);
    observed_co2.T = licor_co2.Tair(:,sample_number);
    observed_co2.E = licor_co2.E(:,sample_number);
    observed_co2.gtc = licor_co2.gtc(:,sample_number);
    observed_co2.Fs = licor_co2.Fs(:,sample_number);
    observed_co2.Fm_ = licor_co2.Fm_(:,sample_number);
    observed_co2.Fm = licor_co2.Fm(:,sample_number);
    observed_co2.Fo = licor_co2.Fo(:,sample_number);
    observed_co2.Fo_ = licor_co2.Fo_(:,sample_number);

    observed_co2.An = licor_co2.A(:,sample_number)*1e-6;
    observed_co2.PhiP = 1- observed_co2.Fs./observed_co2.Fm_;
    observed_co2.PhiN = observed_co2.Fs.*(1./observed_co2.Fm_-1./observed_co2.Fm);
    observed_co2.PhiDF = observed_co2.Fs./observed_co2.Fm; %Y(NO)
    observed_co2.abs = table2array(abs_lrc(sample_number,3));
    observed_co2.ETR = (observed_co2.Q.*0.85./2.*observed_co2.PhiP).*1e-6;

    observed_co2.fs_forw = spec_co2.fs_forw(:,:,sample_number)*7599.14*1e1*2;
    observed_co2.fs_back = spec_co2.fs_back(:,:,sample_number)*7599.14*1e1*2;
    observed_co2.refl = squeeze(spec_co2.refl(:,:,sample_number));
    observed_co2.tran = squeeze(spec_co2.trans(:,:,sample_number));

    % Spectral data
    observed_co2.wl_back_f = spec_co2.wl_back_f;
    observed_co2.wl_forw_f = spec_co2.wl_forw_f;

    observed_co2.smoothed_refl_co2 = zeros(size(observed_co2.refl,1),size(observed_co2.refl,2));
    observed_co2.smoothed_tran_co2 = zeros(size(observed_co2.tran,1),size(observed_co2.tran,2));
    observed_co2.smoothed_Fu_co2 = zeros(size(observed_co2.fs_back,1),size(observed_co2.fs_back,2));
    observed_co2.smoothed_Fd_co2 = zeros(size(observed_co2.fs_forw,1),size(observed_co2.fs_forw,2));
    for n_co2 = 1:15
        observed_co2.smoothed_refl_co2(n_co2,:) = sgolayfilt(observed_co2.refl(n_co2,:),2,9);
        observed_co2.smoothed_tran_co2(n_co2,:) = sgolayfilt(observed_co2.tran(n_co2,:),2,9);
        observed_co2.smoothed_Fu_co2(n_co2,:) =  sgolayfilt(observed_co2.fs_back(n_co2,:),2,9);
        observed_co2.smoothed_Fd_co2(n_co2,:) =  sgolayfilt(observed_co2.fs_forw(n_co2,:),2,9);
    end

    observed_co2.wl_refl_green = spec_co2.wl_refl_green;
    observed_co2.wl_tran_green = spec_co2.wl_trans_green;

    wl_refl_490 = find(observed_co2.wl_refl_green>=490 & observed_co2.wl_refl_green<491);
    wl_refl_520 = find(observed_co2.wl_refl_green>=520 & observed_co2.wl_refl_green<521);
    wl_refl_531 = find(observed_co2.wl_refl_green>=531 & observed_co2.wl_refl_green<532);
    wl_refl_545 = find(observed_co2.wl_tran_green>=545 & observed_co2.wl_refl_green<546);
    wl_refl_570 = find(observed_co2.wl_refl_green>=570 & observed_co2.wl_refl_green<571);

    wl_tran_490 = find(observed_co2.wl_tran_green>=490 & observed_co2.wl_tran_green<491);
    wl_tran_520 = find(observed_co2.wl_tran_green>=520 & observed_co2.wl_tran_green<521);
    wl_tran_531 = find(observed_co2.wl_tran_green>=531 & observed_co2.wl_tran_green<532);
    wl_tran_545 = find(observed_co2.wl_tran_green>=545 & observed_co2.wl_tran_green<546);
    wl_tran_570 = find(observed_co2.wl_tran_green>=570 & observed_co2.wl_tran_green<571);

    refl_490 = observed_co2.smoothed_refl_co2(:,wl_refl_490,:);
    refl_520 = observed_co2.smoothed_refl_co2(:,wl_refl_520,:);
    refl_531 = observed_co2.smoothed_refl_co2(:,wl_refl_531,:);
    refl_545 = observed_co2.smoothed_refl_co2(:,wl_refl_545,:);
    refl_570 = observed_co2.smoothed_refl_co2(:,wl_refl_570,:);

    tran_490 = observed_co2.smoothed_tran_co2(:,wl_tran_490,:);
    tran_520 = observed_co2.smoothed_tran_co2(:,wl_tran_520,:);
    tran_531 = observed_co2.smoothed_tran_co2(:,wl_tran_531,:);
    tran_545 = observed_co2.smoothed_tran_co2(:,wl_tran_545,:);
    tran_570 = observed_co2.smoothed_tran_co2(:,wl_tran_570,:);

    observed_co2.pri_refl = (refl_531-refl_570)./(refl_531+refl_570);
    observed_co2.pri_tran = (tran_531-tran_570)./(tran_531+tran_570);
    observed_co2.tri_pri_refl = 0.5.*((520-490).*(refl_545 - refl_490)-(545-490).*...
        (refl_520-refl_490));
    observed_co2.tri_pri_tran = 0.5.*((520-490).*(tran_545 - tran_490)-(545-490).*...
        (tran_520-tran_490));

    % ---------------------------------------------------------------------
    licor_lrc = observed_real.aus_data.licor_lrc;
    spec_lrc = observed_real.aus_data.spec_lrc;
    observed_lrc.Ca = licor_lrc.Ca(:,sample_number);
    observed_lrc.Ci = licor_lrc.Ci(:,sample_number);
    observed_lrc.Q = licor_lrc.Qin(:,sample_number);
    observed_lrc.T = licor_lrc.Tair(:,sample_number);
    observed_lrc.E = licor_lrc.E(:,sample_number);
    observed_lrc.gtc = licor_lrc.gtc(:,sample_number);
    observed_lrc.Fs = licor_lrc.Fs(:,sample_number);
    observed_lrc.Fm_ = licor_lrc.Fm_(:,sample_number);
    observed_lrc.Fm = licor_lrc.Fm(:,sample_number);
    observed_lrc.Fo = licor_lrc.Fo(:,sample_number);
    observed_lrc.Fo_ = licor_lrc.Fo_(:,sample_number);

    observed_lrc.Q(1) = 5;  % To avoid negative pam
    observed_lrc.Fs(1) = observed_lrc.Fo(1);
    observed_lrc.Fm_(1) = observed_lrc.Fm(1);
    observed_lrc.Fo_(1) = observed_lrc.Fo(1);

    observed_lrc.An = licor_lrc.A(:,sample_number)*1e-6;
    observed_lrc.PhiP = 1- observed_lrc.Fs./observed_lrc.Fm_;
    observed_lrc.PhiN = observed_lrc.Fs.*(1./observed_lrc.Fm_-1./observed_lrc.Fm);
    observed_lrc.PhiDF = observed_lrc.Fs./observed_lrc.Fm; %Y(NO)
    observed_lrc.abs = table2array(abs_lrc(sample_number,3));
    observed_lrc.ETR = (observed_lrc.Q.*0.85./2.*observed_lrc.PhiP).*1e-6;

    observed_lrc.fs_forw = spec_lrc.fs_forw(:,:,sample_number)*7599.14*1e1*2;
    observed_lrc.fs_back = spec_lrc.fs_back(:,:,sample_number)*7599.14*1e1*2;
    observed_lrc.fmax_forw = spec_lrc.fmax_forw(:,:,sample_number)*7599.14*1e1*2;
    observed_lrc.fmax_back = spec_lrc.fmax_back(:,:,sample_number)*7599.14*1e1*2;
    observed_lrc.fmin_forw = spec_lrc.fmin_forw(:,:,sample_number)*7599.14*1e1*2;
    observed_lrc.fmin_back = spec_lrc.fmin_back(:,:,sample_number)*7599.14*1e1*2;

    
    observed_lrc.refl = squeeze(spec_lrc.refl(:,:,sample_number));
    observed_lrc.tran = squeeze(spec_lrc.trans(:,:,sample_number));

    observed_lrc.wl_back_f = spec_lrc.wl_back_f;
    observed_lrc.wl_forw_f = spec_lrc.wl_forw_f;

    observed_lrc.smoothed_refl_lrc = zeros(size(observed_lrc.refl,1),size(observed_lrc.refl,2));
    observed_lrc.smoothed_tran_lrc = zeros(size(observed_lrc.tran,1),size(observed_lrc.tran,2));
    observed_lrc.smoothed_Fu_lrc = zeros(size(observed_lrc.fs_back,1),size(observed_lrc.fs_back,2));
    observed_lrc.smoothed_Fd_lrc = zeros(size(observed_lrc.fs_forw,1),size(observed_lrc.fs_forw,2));
    for n_co2 = 1:16
        observed_lrc.smoothed_refl_lrc(n_co2,:) = sgolayfilt(observed_lrc.refl(n_co2,:),2,9);
        observed_lrc.smoothed_tran_lrc(n_co2,:) = sgolayfilt(observed_lrc.tran(n_co2,:),2,9);
        
        observed_lrc.smoothed_Fu_lrc(n_co2,:) =  sgolayfilt(observed_lrc.fs_back(n_co2,:),2,9);
        observed_lrc.smoothed_Fd_lrc(n_co2,:) =  sgolayfilt(observed_lrc.fs_forw(n_co2,:),2,9);

        observed_lrc.smoothed_Fmax_u_lrc(n_co2,:) =  sgolayfilt(observed_lrc.fmax_back(n_co2,:),2,9);
        observed_lrc.smoothed_Fmax_d_lrc(n_co2,:) =  sgolayfilt(observed_lrc.fmax_forw(n_co2,:),2,9);

        observed_lrc.smoothed_Fmin_u_lrc(n_co2,:) =  sgolayfilt(observed_lrc.fmin_back(n_co2,:),2,9);
        observed_lrc.smoothed_Fmin_d_lrc(n_co2,:) =  sgolayfilt(observed_lrc.fmin_forw(n_co2,:),2,9);
       
    end

    observed_lrc.wl_refl_green = spec_lrc.wl_refl_green;
    observed_lrc.wl_tran_green = spec_lrc.wl_trans_green;


    wl_refl_490 = find(observed_lrc.wl_refl_green>=490 & observed_lrc.wl_refl_green<491);
    wl_refl_520 = find(observed_lrc.wl_refl_green>=520 & observed_lrc.wl_refl_green<521);
    wl_refl_531 = find(observed_lrc.wl_refl_green>=531 & observed_lrc.wl_refl_green<532);
    wl_refl_545 = find(observed_lrc.wl_tran_green>=545 & observed_lrc.wl_refl_green<546);
    wl_refl_570 = find(observed_lrc.wl_refl_green>=570 & observed_lrc.wl_refl_green<571);

    wl_tran_490 = find(observed_lrc.wl_tran_green>=490 & observed_lrc.wl_tran_green<491);
    wl_tran_520 = find(observed_lrc.wl_tran_green>=520 & observed_lrc.wl_tran_green<521);
    wl_tran_531 = find(observed_lrc.wl_tran_green>=531 & observed_lrc.wl_tran_green<532);
    wl_tran_545 = find(observed_lrc.wl_tran_green>=545 & observed_lrc.wl_tran_green<546);
    wl_tran_570 = find(observed_lrc.wl_tran_green>=570 & observed_lrc.wl_tran_green<571);

    refl_490 = observed_lrc.smoothed_refl_lrc(:,wl_refl_490,:);
    refl_520 = observed_lrc.smoothed_refl_lrc(:,wl_refl_520,:);
    refl_531 = observed_lrc.smoothed_refl_lrc(:,wl_refl_531,:);
    refl_545 = observed_lrc.smoothed_refl_lrc(:,wl_refl_545,:);
    refl_570 = observed_lrc.smoothed_refl_lrc(:,wl_refl_570,:);

    tran_490 = observed_lrc.smoothed_tran_lrc(:,wl_tran_490,:);
    tran_520 = observed_lrc.smoothed_tran_lrc(:,wl_tran_520,:);
    tran_531 = observed_lrc.smoothed_tran_lrc(:,wl_tran_531,:);
    tran_545 = observed_lrc.smoothed_tran_lrc(:,wl_tran_545,:);
    tran_570 = observed_lrc.smoothed_tran_lrc(:,wl_tran_570,:);

    observed_lrc.pri_refl = (refl_531-refl_570)./(refl_531+refl_570);
    observed_lrc.pri_tran = (tran_531-tran_570)./(tran_531+tran_570);
    observed_lrc.tri_pri_refl = 0.5.*((520-490).*(refl_545 - refl_490)-(545-490).*...
        (refl_520-refl_490));
    observed_lrc.tri_pri_tran = 0.5.*((520-490).*(tran_545 - tran_490)-(545-490).*...
        (tran_520-tran_490));

    observed.observed_co2 = observed_co2;
    observed.observed_lrc = observed_lrc;
elseif strcmp(dataset,'vilfan')==1
    % Prepare CO2 experiment
    leaf_id_co2 = observed_real.vilfan_data.data_co2.leaf_id;
    leaf_IDs = unique(leaf_id_co2);
    licor_co2 = observed_real.vilfan_data.data_co2.licor_pam;
    I = leaf_id_co2==sample_number;

    observed_co2.Ca = licor_co2.CO2(I);
    observed_co2.Ci = licor_co2.Ci(I);
    observed_co2.PARi = licor_co2.PARi(I);
    observed_co2.T = licor_co2.Tleaf(I);
    observed_co2.RHS = licor_co2.RHS(I);
    observed_co2.Fs = licor_co2.Ft(I);
    observed_co2.Fm_ = licor_co2.Fmp(I);
    observed_co2.Fm = licor_co2.Fm(I);
    observed_co2.Fo = licor_co2.Fo(I);
    observed_co2.PRI = licor_co2.PRI(I);
    observed_co2.Y = licor_co2.Y(I);
    observed_co2.Vcmax = licor_co2.Vcmax(I);

    observed_co2.An = licor_co2.Photo(I)*1e-6;
    observed_co2.PhiP = 1- observed_co2.Fs./observed_co2.Fm_;
    observed_co2.PhiN = observed_co2.Fs.*(1./observed_co2.Fm_-1./observed_co2.Fm);
    observed_co2.PhiDF = observed_co2.Fs./observed_co2.Fm;
    %     observed_co2.abs = table2array(abs_lrc(sample_number,3));
    observed_co2.ETR = (observed_co2.PARi.*0.85./2.*observed_co2.PhiP).*1e-6;

    observed_co2.fs_forw = observed_real.vilfan_data.data_co2.Fd(:,I');
    observed_co2.If = observed_real.vilfan_data.data_co2.If(:,I');
    observed_co2.Ti = observed_real.vilfan_data.data_co2.Ti(:,I');
    observed_co2.Tnorm565 = observed_real.vilfan_data.data_co2.TnormRatio565(:,I');
    observed_co2.wl = (350:1:2500)';

    observed_co2.smoothed_fs_forw = zeros(size(observed_co2.fs_forw,1),size(observed_co2.fs_forw,2));
    observed_co2.smoothed_If = zeros(size(observed_co2.If,1),size(observed_co2.If,2));
    observed_co2.smoothed_Ti = zeros(size(observed_co2.Ti,1),size(observed_co2.Ti,2));
    observed_co2.smoothed_Tnorm565 = zeros(size(observed_co2.Tnorm565,1),size(observed_co2.Tnorm565,2));
    for n_co2 = 1:10
        observed_co2.smoothed_fs_forw(:,n_co2) = sgolayfilt(observed_co2.fs_forw(:,n_co2),2,9);
        observed_co2.smoothed_If(:,n_co2) = sgolayfilt(observed_co2.If(:,n_co2),2,9);
        observed_co2.smoothed_Ti(:,n_co2) = sgolayfilt(observed_co2.Ti(:,n_co2),2,9);
        observed_co2.smoothed_Tnorm565(:,n_co2) = sgolayfilt(observed_co2.Tnorm565(:,n_co2),2,9);

    end


    % Prepare LRC experiment
    leaf_id_lrc = observed_real.vilfan_data.data_lrc.leaf_id;
    leaf_IDs = unique(leaf_id_lrc);
    licor_lrc = observed_real.vilfan_data.data_lrc.licor_pam;
    I = leaf_id_lrc==sample_number;

    observed_lrc.Ca = licor_lrc.CO2(I);
    observed_lrc.Ci = licor_lrc.Ci(I);
    observed_lrc.PARi = licor_lrc.PARi(I);
    observed_lrc.T = licor_lrc.Tleaf(I);
    observed_lrc.RHS = licor_lrc.RHS(I);
    observed_lrc.Fs = licor_lrc.Ft(I);
    observed_lrc.Fm_ = licor_lrc.Fmp(I);
    observed_lrc.Fm = licor_lrc.Fm(I);
    %     observed_lrc.Fo = licor_lrc.Fo(I);
    observed_lrc.PRI = licor_lrc.PRI(I);

    observed_lrc.An = licor_lrc.Photo(I)*1e-6;
    observed_lrc.PhiP = 1- observed_lrc.Fs./observed_lrc.Fm_;
    observed_lrc.PhiN = observed_lrc.Fs.*(1./observed_lrc.Fm_-1./observed_lrc.Fm);
    observed_lrc.PhiDF = observed_lrc.Fs./observed_lrc.Fm;
    %     observed_lrc.abs = table2array(abs_lrc(sample_number,3));
    observed_lrc.ETR = (observed_lrc.PARi.*0.85./2.*observed_lrc.PhiP).*1e-6;

    observed_lrc.fs_forw = observed_real.vilfan_data.data_lrc.Fd(:,I');
    observed_lrc.If = observed_real.vilfan_data.data_lrc.If(:,I');
    observed_lrc.Ti = observed_real.vilfan_data.data_lrc.Ti(:,I');
    observed_lrc.Tnorm565 = observed_real.vilfan_data.data_lrc.TnormRatio565(:,I');
    observed_lrc.wl = (350:1:2500)';

    observed_lrc.smoothed_fs_forw = zeros(size(observed_lrc.fs_forw,1),size(observed_lrc.fs_forw,2));
    observed_lrc.smoothed_If = zeros(size(observed_lrc.If,1),size(observed_lrc.If,2));
    observed_lrc.smoothed_Ti = zeros(size(observed_lrc.Ti,1),size(observed_lrc.Ti,2));
    observed_lrc.smoothed_Tnorm565 = zeros(size(observed_lrc.Tnorm565,1),size(observed_lrc.Tnorm565,2));
    for n_lrc = 1:9
        observed_lrc.smoothed_fs_forw(:,n_lrc) = sgolayfilt(observed_lrc.fs_forw(:,n_lrc),2,9);
        observed_lrc.smoothed_If(:,n_lrc) = sgolayfilt(observed_lrc.If(:,n_lrc),2,9);
        observed_lrc.smoothed_Ti(:,n_lrc) = sgolayfilt(observed_lrc.Ti(:,n_lrc),2,9);
        observed_lrc.smoothed_Tnorm565(:,n_lrc) = sgolayfilt(observed_lrc.Tnorm565(:,n_lrc),2,9);

    end
    observed.observed_co2 = observed_co2;
    observed.observed_lrc = observed_lrc;
end