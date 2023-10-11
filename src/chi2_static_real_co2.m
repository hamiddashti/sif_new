function chi2=chi2_static_real_co2(vars,obs,weights,v_inv)
% v_inv.e1e2 =vars(1); 
v_inv.Vqmax = vars(1);
v_inv.Vcmax = vars(2);
v_inv.rm = vars(3);
% v_inv.eps1 = v_inv.e1e2./(1+v_inv.e1e2);
% v_inv.eps2 = 1./(1+v_inv.e1e2);


an_obs = obs.An;
PhiP_obs = obs.PhiP;
PhiN_obs = obs.PhiN;
PhiDF_obs = obs.PhiDF;
% ETR_obs = obs.ETR;

wAn = weights.wAn;
wPhiP = weights.wPhiP;
wPhiN = weights.wPhiN;
wPhiDF = weights.wPhiDF;
% wETR = weights.wETR;

try
m_syn=model_fun_static(v_inv);
an_sim = m_syn.An_a;
PhiP_sim = m_syn.PAM1_a;
PhiN_sim = m_syn.PAM2_a;
PhiDF_sim = m_syn.PAM3_a;
% ETR_sim = m_syn.PAM4_a;

an_obs_norm = (an_obs-mean(an_obs))./std(an_obs);
an_sim_norm = (an_sim-mean(an_obs))./std(an_obs);

PhiP_obs_norm = (PhiP_obs-mean(PhiP_obs))./std(PhiP_obs);
PhiP_sim_norm = (PhiP_sim-mean(PhiP_obs))./std(PhiP_obs);

PhiN_obs_norm = (PhiN_obs-mean(PhiN_obs))./std(PhiN_obs);
PhiN_sim_norm = (PhiN_sim-mean(PhiN_obs))./std(PhiN_obs);

PhiDF_obs_norm = (PhiDF_obs-mean(PhiDF_obs))./std(PhiDF_obs);
PhiDF_sim_norm = (PhiDF_sim-mean(PhiDF_obs))./std(PhiDF_obs);
% 
% ETR_obs_norm = (ETR_obs-mean(ETR_obs))./std(ETR_obs);
% ETR_sim_norm = (ETR_sim-mean(ETR_obs))./std(ETR_obs);


a(1) = mean(abs(an_sim_norm-an_obs_norm));
a(2) = mean(abs(PhiP_sim_norm-PhiP_obs_norm));
a(3) = mean(abs(PhiN_sim_norm-PhiN_obs_norm));
a(4) = mean(abs(PhiDF_sim_norm-PhiDF_obs_norm));
% a(5) = mean(abs(ETR_sim_norm-ETR_obs_norm));

chi2 = (wAn*a(1))+(wPhiP*a(2))+(wPhiN*a(3))+(wPhiDF*a(4));
catch
% disp(v_inv)
chi2 = NaN;
end