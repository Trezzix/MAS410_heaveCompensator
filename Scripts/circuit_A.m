clc;clear;close all;
% format bank

%%%%%%%%%%%%%%%%%%%%% Given Constants %%%%%%%%%%%%%%%%%%%%%

g = 9.81;  % [m/s^2]
dD = 0.45; % [m], diameter Drum
dR = 0.5;  % [m], diameter gear rim
dp = 0.15; % [m], diameter pinion
mu_eq = 0.15; % equivalent friction coefficient
w0 = 5; % [rad/s]
ig = 7; % gear ratio motor -> pinion
n_sh = 3; % number of sheaves
mpl = 24000; % [kg], payload mass
Zw = 1.2;  % [m], wave amplitude
Tw = 10.0; % [s], wave period
eta_hmM = 1;   % hydromechanical effiency of motor
eta_vM = 0.94; % volumetric effiency of motor

%%%%%%%%%%%%%%%%%%%%% Chosen Constants %%%%%%%%%%%%%%%%%%%%%

ps = 230e5; % 210 [bar] -> [Pa]
nm = 2; % number of motors
    % valve
nv = 1; % number of valves
pr = 10e5; % [bar] -> [Pa] ,datasheet
safety_factor = 1.1; % from lecture
    % TODO: Liquid Properties
rho = 875; % [kg/m^3] placeholder Liquid Density
beta = 900e6; % [MPa] -> [Pa] placeholder Liquid Stiffness (Bulk Modulus)

%%%%%%%%%%%%%%%%% General for all configs %%%%%%%%%%%%%%%%%
% gear ratios
i_p = (dR/2)/(dp/2);
iT = ig * i_p; % total gear ratio between motor and drum
i_pl2M = (dD/2)/(2*n_sh*ig*i_p*nm); % [m]

% Max Speed
    % don't use transmission from payload2motor, since all motors
    % retain the same speed (thus don't multiply with nm)
zdot_max = (Zw * 2*pi) / Tw; % t=Tw --> cos(2pi)=1
thetadot_D_max = zdot_max*(2*n_sh*2/dD); % max speed of drum [rad/s], absolute value, wrong?
thetadot_m_max = iT*thetadot_D_max; % [rad/s]
thetadot_m_max_rpm = thetadot_m_max * (60/(2*pi)); % [rad/s] -> [RPM]
MaxSpeed = table(thetadot_D_max, thetadot_m_max, thetadot_m_max_rpm)

% Max Acceleration
zdotdot_max = -((Zw*(2*pi)^2)/(Tw^2)); % [m/s^2]
thetadotdot_D_max = -(4*n_sh*zdotdot_max)/dD; % [rad/s^2]
thetadotdot_m_max = iT*thetadotdot_D_max; % [rad/s^2]
MaxAccel = table(thetadotdot_D_max, thetadotdot_m_max)

%%%%%%%%%%%%%%%%%%%%% Chosen Specific %%%%%%%%%%%%%%%%%%%%%
% Max Motor Moment
M_M_max = ((mpl * g * dD * dp) / (4 * n_sh * dR * ig * nm)) * ...
          (1 + mu_eq * tanh(thetadot_m_max/w0))

% Largest motor we can use with speed < calculated above (general for all)
D_min_rpm = 90; % [RPM]

% For choosing motor size
pL_assume = (2/3) * ps; % chosen
D_min = (2*pi * M_M_max) / pL_assume; % ps -> delta_m_p(?) 
    % TODO: might be missing volumetric efficiency
D_min_cm = D_min * 1e6;
MinDisplacement = table(D_min, D_min_cm, D_min_rpm)
% Auto Choosing Motor Size
motorType = [4.93 10.3 12 16 22.9 28.1 32 45.6 56.1 63 80.4 ...
             90 106.7 125 160.4 180 200 250 355 500 710 1000];
MotorJ = [0.00006 0.0004 0.0004 0.0004 0.0012 0.0012 0.0012 0.0024 0.0042 0.0042 0.0072 0.0072 ...
    0.0116 0.0116 0.0220 0.0220 0.0353 0.061 0.102 0.178 0.55 0.55]; % [kg/m^2]
if D_min_cm > D_min_rpm
    warning("Motor size greater than 90, increase nm")
end

for i_for = 1:length(motorType)
    if motorType(i_for) > D_min_cm
        Dm_cm = motorType(i_for);
        Dm = motorType(i_for) * 1e-6;
        Jm = MotorJ(i_for);
        break
    end
end
Jtot = Jm + (mpl)*i_pl2M^2; % TODO: double check
chosenMotor = table(Dm_cm, Dm, Jm, Jtot)

% pressure
pL_max = (M_M_max + Jtot * thetadotdot_m_max) * ((2*pi)/Dm);
% pL = (M_M_max) * ((2*pi)/Dm);
pL_max_bar = pL_max*1e-5

%%%%%%%%%%%%%%%%%%%%% Circuit Specific %%%%%%%%%%%%%%%%%%%%%
% Theoretical Flow
Qm_t = (Dm/(2*pi)) * thetadot_m_max; % [m^3/sec]
Qm_t_Lpmin = Qm_t * 6*10^4;          % [L/min]
% No-Load Flow
Qm_NL = Qm_t * sqrt(ps/(ps - pL_max)); % [m^3/sec]
Qm_NL_Lpmin = Qm_NL * 6*10^4;          % [L/min]
% Leakage Flow
QL = (Qm_t*(1-eta_vM))/eta_vM;    % [m^3/sec]
QL_Lpmin = QL * 6*10^4;           % [L/min]
CdAd_L = QL/sqrt((2/rho)*pL_max); % [L/min]

% servo valve sizing
Qm_NL_total = ((Qm_NL)*nm)/(nv);
Qr_min = safety_factor * Qm_NL_total * sqrt(pr/ps);
Qr_min_lpmin = Qr_min * 6*10^4 % [L/min];

servoFlowsTheoretical = table(Qm_t_Lpmin, Qm_NL_Lpmin, QL_Lpmin)

% Auto Choosing Valve Size
servoType = [30 60 80 150 250 350 550 1000 1500]; % [L/min]
for i_for = 1:length(servoType)
    if servoType(i_for) > Qr_min_lpmin
        Qr_lpmin = servoType(i_for);
        Qr = servoType(i_for) / (6*10^4); % [L/min] -> [m^3/s]
        break
    end
end
servoFlows = table(Qr_min_lpmin, Qr_lpmin, Qr)

CdAd = Qr/sqrt((2/rho) * pr);
Cd = 0.6;
Ad = CdAd/Cd;
servoCdAd = table(Cd, Ad)