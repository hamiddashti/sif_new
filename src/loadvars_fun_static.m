function [Q, T, O, ...
    Abs, beta, CB6F, RUB, Rds, ...
    Kf, Kd, Kp1, Kn1, Kp2, Ku2, ...
    kq, nl, nc, kc, ko, Kc, Ko, ...
    alpha_opt,eps1, eps2,...
    Ci,ss,rm,Vqmax,Vcmax,a2]= loadvars_fun_static(v)

% Rename default values of parameters passed in from 'v' structure
alpha_opt = v.alpha_opt;
a2 = NaN;
Ci = v.Ci;
ss = v.ss;

Vqmax = v.Vqmax;
Vcmax = v.Vcmax;
rm=v.rm;

% Environmental variables
Q = v.Q;
T = v.T;
O = v.O;

% Physiological variables
Abs = v.Abs;
beta = v.beta;
CB6F = v.CB6F;
RUB = v.RUB;
Rds = v.Rds;

% Photochemical constants
Kf = v.Kf;
Kd = v.Kd;
Kp1 = v.Kp1;
Kn1 = v.Kn1; 
Kp2 = v.Kp2;
Ku2 = v.Ku2;

% Biochemical constants
kq = v.kq;
nl = v.nl;
nc = v.nc;
kc = v.kc;
ko = v.ko;
Kc = v.Kc;
Ko = v.Ko;

% solve_xcs = v.solve_xcs;
theta1 = v.theta1;
eps1 = v.eps1;
eps2 = v.eps2;

end