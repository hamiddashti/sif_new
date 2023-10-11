function dyno_sols = dynamic_solver()

%% Symbolic solutions for C3 model

% Rev 2023-03-10
% Definition of symbols
% Clean up working environment
clearvars -except currdir workdir datadir resultsdir; % variables
warning('off'); % suppress warnings

% Define symbols
syms C_steps O Kc Ko kc ko S Q CB6F RUB kq Kp1 Kp2 Kd Kf Ku nl nc...
    JP700_j JP700_c a2 a1 Abs phi2P_a phi2p_a ...
    phi2u_a q2_a JP680_j Vc_j An_j Ag_j An_c Ag_c E gtc gm Ca Vc_c Rd;

S = (kc./Kc).*(Ko./ko); % Rubisco specificity for CO2/O2, dimensionless

gammas = O./(2.*S); % CO2 compensation point in the absence of Rd, bar

% -------------------------------------------------------------------- %

% (1) If Ca known + STs active, calculate a2 as f(Cc_j):

% -------------------------------------------------------------------- %

eq1 = JP700_j - (Q.*CB6F.*kq)./(Q+(CB6F.*kq)./((Abs-a2).*(Kp1./(Kp1 + Kd + Kf)))) == 0;

eq2 = subs(eq1, JP700_j, JP680_j.*(1-(nl./nc)+(3+7.*gammas./C_steps)./((4+8.*gammas./C_steps).*nc)));

eq3 = subs(eq2, JP680_j, Vc_j.*(4.*(1+2.*gammas./C_steps)));

eq4 = subs(eq3, Vc_j, Ag_j./(1 - gammas./C_steps));

eq5 = subs(eq4, Ag_j, An_j + Rd);

eq6 = subs(eq5, An_j, -(C_steps.*E.*gm + Ca.*E.*gm + 2.*C_steps.*gm.*gtc - 2.*Ca.*gm.*gtc)./(E + 2.*gm + 2.*gtc));

eq7 = solve(eq6, a2, 'MaxDegree',3); % a2 as f(Cc)

clear eq1 eq2 eq3 eq4 eq5 eq6;

% Pass out first and only root:

eqA = matlabFunction(eq7(1,1));

% fn is @(Abs,C,CB6F,Ca,E,Kc,Kd,Kf,Ko,Kp1,O,Q,Rd,gm,gtc,kc,ko,kq,nc,nl)

% -------------------------------------------------------------------- %

% (2) If Ca known + STs active, calculate a2 as f(Cc_j):

% -------------------------------------------------------------------- %

eq8 = JP700_j - Q.*a2.*phi2P_a.*(1-(nl./nc)+(3+7.*gammas./C_steps)./((4+8.*gammas./C_steps).*nc)) == 0;

eq9 = subs(eq8,JP700_j, Q.*CB6F.*kq./(Q+CB6F.*kq./((Abs-a2).*(Kp1./(Kp1 + Kd + Kf)))));

eq10 = subs(eq9, phi2P_a, phi2p_a./(1-phi2u_a));

eq11 = subs(eq10, phi2p_a, (q2_a).*Kp2./(Kp2 + Kd + Kf + Ku));

eq12 = subs(eq11, phi2u_a, (q2_a).*Ku./(Kp2 + Kd + Kf + Ku) + (1 - q2_a).*Ku./(Kd + Kf + Ku));

eq13 = subs(eq12, q2_a, 1 - (Q./(Q+CB6F.*kq./((Abs-a2).*(Kp1./(Kp1 + Kd + Kf))))));

eq14 = solve(eq13, a2, 'MaxDegree',3); % a2 as f(Cc)

clear eq8 eq9 eq10 eq11 eq12 eq13;

% Pass out second root:

eqB = matlabFunction(eq14(2,1));

% fn is @(Abs,C,CB6F,Kc,Kd,Kf,Ko,Kp1,Kp2,Ku,O,Q,kc,ko,kq,nc,nl)

% -------------------------------------------------------------------- %

% (3) Plot both functions numerically and find intersection

% -------------------------------------------------------------------- %

clearvars -except eqA eqB
dyno_sols = workspace2struct_fun();
end
