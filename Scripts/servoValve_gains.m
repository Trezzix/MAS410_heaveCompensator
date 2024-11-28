clc; clear; close all

%% Servo Valve Gains

% From circuit_A / simulink model
ps = 220e5; % [Pa]
VL1 = 1e-3; % [m^3]
VL2 = VL1;
beta = 1000e6; % [Pa]
rho = 875; % [kg/m^3]
Cd = 0.6;
pr = 10e5; % [Pa]
Qr = 0.0092; % [m^3/sec]
QL_ss = 0.0047; % [m^3/sec]
pL_ss = 166.6e5; % [Pa]
Dm = 6.3e-5; % [m^3/rev]
Dw = Dm/(2*pi); % [m^3/rad]
% Jpl = 0.062; % [kg*m^2]
J = 0.074589795918367;

Vtot = VL1 + VL2;
% Hydraulic Capacitance
C = Vtot / (4*beta);
Ad0 = Qr / (Cd * sqrt(pr/rho) ); % [m^2]

u_ss = QL_ss / (Cd*Ad0*sqrt( (ps - pL_ss)/rho ) );

K_qu = Cd*Ad0*sqrt( ( ps - sign(u_ss)*pL_ss ) / rho );
K_qp = Cd*Ad0*abs(u_ss) / (2 * sqrt( rho*(ps - sign(u_ss)*pL_ss) ) );

K_mh = K_qu / Dw;

omega_mh = Dw / sqrt(J*C);
zeta_mh = K_qp / (2*Dw) * sqrt(J/C);

% Assume infinitely fast valve?
omega_vmh = 1*omega_mh;
Kp = 0.7079 * (2*zeta_mh*omega_vmh) / K_mh; % [1/m]