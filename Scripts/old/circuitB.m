clc;clear;close all;


% Chosen Constants
ps = 190e5; % 210 [Pa] Supply Pressure
nM = 1; % Number of Motors

% From datasheet - Size 125
% J_motor = 0.0116; % [kg*m^2]
DM = 125e6; % [m^3/rev]
DM_cm = DM*1e-6; % [cm^3/rev]

% Cost [-]
wM = 2;
DMmax = 1000e-3; % [m^3/rev]
cost_motor = wM * (1 + DM/DMmax);
cost_motors = cost_motor * nM;

table(cost_motor, cost_motors, DM_cm)

motorType = [4.93 10.3 12 16 22.9 28.1 32 45.6 56.1 63 80.4 ...
             90 106.7 125 160.4 180 200 250 355 500 710 1000];
for idx = 1 : length(motorType)
    costMotors(idx) = wM*(1 + motorType(idx)/DMmax);
end

figure; hold on;
for i = 1 : 5
    plot(costMotors*i)
end
legend('1 Motor','2 Motors','3 Motors','4 Motors','5 Motors', ...
       'location','northwest')
title('Cost of Motors')
xlabel('Motor Size #')
ylabel('Cost [-]')

%% 

% Given Constants
g = 9.81; % [m/s^2]
dD = 0.45; % [m] Diameter Drum
dR = 0.5; % [m] Diameter Gear Rim
dp = 0.15; % [m] Diameter Pinion
mu_eq = 0.15; % Equivalent Friction Coefficient
ig = 7; % Gear Ratio Motor -> Pinion
n_sh = 3; % Number of Sheaves
m_pl = 24000; % [kg] Payload Mass
eta_vm = 1; % Hydromechanical Effiency of Motor

w0 = 5; % [rad/sec]
Zw = 1.2; % [m] Wave Amplitude
Tw = 10.0; % [sec] Wave Period

% Calculations
zDot_max = Zw * 2*pi / Tw * 1; % [m/s]

M_M_max = ((m_pl * g * dD * dp) / (4 * n_sh * dR * ig * nM)) * ...
          (1 + mu_eq * tanh(thetadot_m_max/w0)); % When hoisting

D_min = (2*pi * M_M_max) / ps; %ps -> delta_m_p(?)
D_min_cm = D_min * 1e6
thetadot_m_max_rpm = thetadot_m_max * (2*pi)
%T = (D_min_cm * ps *1e-5 * eta_vm)/ (20 * pi) %singular motor, from datasheet, uses delta_p
%flow
Qm_max = D_min * thetadot_m_max; %WIP: D NEEDS TO BE CONVERTED TO Dw
Qm_max_L_min = Qm_max * 6*10^4


