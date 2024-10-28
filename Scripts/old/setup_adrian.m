clc;clear;close all;
%given constants
g = 9.81; % m/s^2
dD = 0.45; %m, diameter Drum
dR = 0.5; %m, diameter gear rim
dp = 0.15; %m, diameter pinion
mu_eq = 0.15; %equivalent friction coef
w0 = 5; %rad/s, 
ig = 7; % gear ratio motor -> pinion
n_sh = 3; %number of sheaves
mpl = 24000; %kg, payload mass
Zw = 1.2; %m, wave amplitude
Tw = 10.0; %s, wave period
eta_vm = 1; %hydromechanical effiency of motor

%chosen constants
ps = 190e5; %210 bar -> Pa
nm = 1; %number of motors

%calculating max speed
zdot_max = (Zw * 2*pi) / Tw; %absolute value
thetadot_D_max = (12*zdot_max)/dD %max speed of drum [rad/s]
thetadot_m_max = ((dR/2)*ig*thetadot_D_max)/(dp/2) %[rad/s]
M_M_max = ((mpl * g * dD * dp) / (4 * n_sh * dR * ig * nm)) * (1 + mu_eq * tanh(thetadot_m_max/w0))
D_min = (2*pi * M_M_max) / ps; %ps -> delta_m_p(?)
D_min_cm = D_min * 1e6
thetadot_m_max_rpm = thetadot_m_max * (2*pi)
%T = (D_min_cm * ps *1e-5 * eta_vm)/ (20 * pi) %singular motor, from datasheet, uses delta_p
%flow
Qm_max = D_min * thetadot_m_max; %WIP: D NEEDS TO BE CONVERTED TO Dw
Qm_max_L_min = Qm_max * 6*10^4