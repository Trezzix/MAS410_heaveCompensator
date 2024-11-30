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
ps = 220e5; % [bar] -> [Pa]
nm = 3; % number of motors
    % Proportional Valve
npv = 1; % number of proportional valves
% pr = 10e5; % [bar] -> [Pa], update spool flows below from datasheet
pN = 350e5; % [Pa] Nominal pressure from datasheet
deltaP_spool = 4.2e5; % [bar] -> [Pa], from datasheet of cvg50 31-08
deltaP_comp = 6e5; % [bar] -> [Pa], pg 12, cvg50
    % Counterbalance Valve (CBV)
ncbv = 2; % number of CBV
pB = 18e5; % choose between 10-30 [bar], higher = bigger pressure in p1, lower = cavitation
pcr2_over = 1.3; % 1.1 to 1.3 (see tutorial 4)
    % Liquid Properties
rho = 875; % [kg/m^3] Liquid Density
beta = 1000e6; % [MPa] -> [Pa] Liquid Stiffness (Bulk Modulus)

% Failure variables
eta_vM_failure = 0.50; % 50 percent

% [WARNING] Remember to update datasheet Qr from lines 130 onwards!

%%%%%%%%%%%%%%%%% General for all configs %%%%%%%%%%%%%%%%%
% gear ratios
i_p = (dR/2)/(dp/2);
iT = ig * i_p; % total gear ratio between motor and drum
i_pl2M = (dD/2)/(2*n_sh*ig*i_p); % [m]

% Max Speed
    % don't use transmission from payload2motor, since all motors
    % retain the same speed (thus don't multiply with nm)
zdot_max = (Zw * 2*pi) / Tw; % t=Tw --> cos(2pi)=1
thetadot_D_max = zdot_max*(2*n_sh*2/dD); % max speed of drum [rad/s]
thetadot_m_max = iT*thetadot_D_max; % [rad/s]
thetadot_m_max_rpm = thetadot_m_max * (60/(2*pi)); % [rad/s] -> [RPM]
MaxSpeed = table(thetadot_D_max, thetadot_m_max, thetadot_m_max_rpm)

% Max Acceleration
zdotdot_max = -((Zw*(2*pi)^2)/(Tw^2)); % [m/s^2]
thetadotdot_D_max = -(4*n_sh*zdotdot_max)/dD; % [rad/s^2]
thetadotdot_m_max = iT*thetadotdot_D_max; % [rad/s^2]
MaxAccel = table(thetadotdot_D_max, thetadotdot_m_max)

%%%%%%%%%%%%%%%%%%%%% Chosen Specific %%%%%%%%%%%%%%%%%%%%%
    % Volumetric Efficiency taken into account in Leakage flow below
    % Hydromechanical efficiency is taken into account, project eq. (2)

% Max Motor Moment
M_M_max = ((mpl * g * dD * dp) / (4 * n_sh * dR * ig * nm)) * ...
          (1 + mu_eq * tanh(thetadot_m_max/w0))

% Largest motor we can use with speed < calculated above (general for all)
Dmax_minRPM = 90; % [cm^3/rev]

% For choosing motor size
pL_assume = (2/3) * ps; % chosen, different between circuits
D_min = (2*pi * M_M_max) / pL_assume;
D_min_cm = D_min * 1e6;
MotorDisplacements = table(D_min, D_min_cm, Dmax_minRPM)
% Auto Choosing Motor Size
motorType = [4.93 10.3 12 16 22.9 28.1 32 45.6 56.1 63 80.4 ...
             90 106.7 125 160.4 180 200 250 355 500 710 1000];
MotorJ = [0.00006 0.0004 0.0004 0.0004 0.0012 0.0012 0.0012 ...
          0.0024 0.0042 0.0042 0.0072 0.0072 0.0116 0.0116 ...
          0.0220 0.0220 0.0353 0.061 0.102 0.178 0.55 0.55]; % [kg/m^2]
if D_min_cm > Dmax_minRPM
    error("Motor size greater than 90, increase nm")
end

% Find smallest motor above minimum
for i_for = 1:length(motorType)
    if motorType(i_for) > D_min_cm
        Dm_cm = motorType(i_for);
        Dm = motorType(i_for) * 1e-6;
        Jm = MotorJ(i_for);
        break
    end
end
Jpl = (mpl)*i_pl2M^2; % [kg*m^2], payload inertia
Jtot = Jm*nm + Jpl;      % [kg*m^2], total inertia
chosenMotor = table(Dm_cm, Dm, Jm, Jpl, Jtot)

% pressure
% pL_max = (M_M_max + Jtot * thetadotdot_m_max) * ((2*pi)/Dm);
pL_max = (M_M_max) * ((2*pi)/Dm);
pL_max_bar = pL_max*1e-5

%%%%%%%%%%%%%%%%%%%%% Circuit B Specific %%%%%%%%%%%%%%%%%%%%%
% Theoretical Flow
Qm_t = (Dm/(2*pi)) * thetadot_m_max; % [m^3/sec]
Qm_t_Lpmin = Qm_t * 6*10^4; % [L/min], for proportional valve sizing
% Leakage Flow
pL_leakage_motor_max = 350e5; % [bar] from motor datasheet
QL = (Qm_t*(1-eta_vM))/eta_vM;    % [m^3/sec]
QL_Lpmin = QL * 6*10^4;           % [L/min]
CdAd_L = QL/sqrt((2/rho)*pL_leakage_motor_max); % [m^2]
Cd_L = 0.6;
Ad_L = CdAd_L/Cd_L; % [m^2]
% Failure - High Leakage Flow
QL_fail = (Qm_t*(1-eta_vM_failure))/eta_vM_failure;    % [m^3/sec]
% QL_Lpmin = QL * 6*10^4;           % [L/min]
CdAd_L_fault = QL_fail/sqrt((2/rho)*pL_leakage_motor_max); % [m^2]
Ad_L_fault = CdAd_L_fault/Cd_L; % [m^2]

leakageFlows = table(Cd_L, Ad_L, Ad_L_fault)

%%%%% Proportional Valve - Spool %%%%%
% type: closed center, symmetric, CVGxx 31-xx, datasheet page 16
Qm_max_total = ((Qm_t+QL)*nm)/(npv);
Qm_max_total_lpmin = Qm_max_total * 6*10^4;
spoolFlows = [60, 140, 220, 500, 820, 1000, 950, 1150]; % [L/min]
spoolTypes = ["CVG30 31-00", "CVG30 31-01", "CVG30 31-02", "CVG30 31-05",...
    "CVG50 31-08", "CVG50 31-10", "CVG60 31-10", "CVG60 31-20"];

% Valve must supply at least worst case scenario (max Q from all motors)
for i_for = 1:length(spoolFlows)
    if spoolFlows(i_for) > Qm_max_total_lpmin
        Q_nom_spool = spoolFlows(i_for);
        chosenSpool = spoolTypes(i_for);
        break
    end
end
table(Q_nom_spool, chosenSpool)

% Datasheet: Read pressure @ Qm_max_total [L/min]
prMain = 3.3e5; % [Pa] from datasheet - CVG50 31-10 page 14
prComp = 6.5e5; % [Pa] from datasheet - CVG50 31-10 page 12

% Coefficients and areas for main & compensator spools
CdAd_mainSpool = Qm_max_total/sqrt((2/rho) * prMain); % [m^2]
Cd = 0.6;
Ad = CdAd_mainSpool/Cd;
CdAd_compSpool = Qm_max_total/sqrt((2/rho) * prComp); % [m^2]
Cd_comp = Cd;
Ad_comp = CdAd_compSpool/Cd_comp;

% Compensator Spring
% pcr1 = (Qm_max_total^2 * rho) / (CdAd_mainSpool^2 * 2); % 3.36
% pcr1_bar = pcr1 * 1e-5;
% 
% table(CdAd_mainSpool, CdAd_compSpool, pcr1_bar)

% delta_p_spool = (Qm_max_total^2 * rho) / (CdAd_spool^2 * 2);
% pM_in = ps - (deltaP_comp + deltaP_spool); % [bar], into motor

%%%%% Counterbalance / Overbalance Valve %%%%%
    % Only during lowering, otherwise only goes through check valve
cbv_alpha_list = [1.5 2 2.3 3 4.5 10];
cbv_name_list = ["CBIB" "CBIY" "CBIL" "CBIA" "CBIG" "CBIH"];
% pM_out = pM_in - pL_max; % unsure if correct, TODO: double check
% CBV Spring
pcr2 = pL_max * pcr2_over; % pL_max is very close to the method used in the examples
pcr2_bar = pcr2 * 1e-5;
% Return pressure: A-side Valve
pTank = 0e5; % [Pa]
pA = (Qm_max_total^2 * rho)/(CdAd_mainSpool^2 * 2) + pTank; % 3.36
pA_bar = pA * 1e-5;
% Motor A-side Pressure
pAm = ((M_M_max * 2 * pi) / Dm) + pB;
    % Pilot Ratio
alpha_max = (pAm - pcr2 - pA) / (pA - pB);
% alpha_max = (pL_max + pB - pcr2 - pA) / (pA - pB);

% Smallest CBV with sufficient alpha
% for i_for = 1:length(cbv_alpha_list)
%     if cbv_alpha_list(i_for) > alpha_max
%         alpha_cbv = cbv_alpha_list(i_for-1);
%         cbv_type = cbv_name_list(i_for-1);
%         break
%     end
% enda
diffList_CBValpha = cbv_alpha_list - alpha_max;
[~, alphaIDX] = min(abs(diffList_CBValpha));
alpha_cbv = cbv_alpha_list(alphaIDX);
cbv_type = cbv_name_list(alphaIDX);

cbvStats = table(pcr2_bar, pA_bar, alpha_max, alpha_cbv, cbv_type)

% pB = (pAm - pcr2 + pA*(-1 -alpha_cbv))/((-1 -alpha_cbv));
pB = (pcr2 + (1+alpha_cbv)*pA - pAm);
pB_bar = pB * 1e-5
max_capacity = 480; % [L/min]
n_cbv_min = ceil(Qm_max_total_lpmin/max_capacity)

% Q_free_flow_chk = 320 / 6e4; % [l/min]
% Q_free_flow_cbv = 330 / 6e4; % [l/min]
Q_free_flow_chk = 480 / 6e4; % [l/min]
Q_free_flow_cbv = 480 / 6e4; % [l/min]
pCHK = 22e5; % [Pa]
pCBV = 30e5; % [Pa]
CdAd_chk =      Q_free_flow_chk/sqrt((2/rho) * pCHK); % [m^2]
CdCHK = Cd;
AdCHK = CdAd_chk/CdCHK;
CdAd_cbv_free = Q_free_flow_cbv/sqrt((2/rho) * pCBV); % [m^2]
CdCBV = Cd;
AdCBV = CdAd_cbv_free/CdCBV;

%%%%%%%%%%%%%%%%%%%%%%%% Valve Dynamics %%%%%%%%%%%%%%%%%%%%%%%%
% Read datasheet CVG50 step response - page 12

Tp = 390e-3; % [sec]
cMax   = 14.0001; % [mm]
cFinal = 14.0; % [mm]
OS = (cMax - cFinal)/cFinal;
zeta = (-log(OS))/sqrt(pi^2 + (log(OS))^2); % log = ln for MATLAB :)
wn = pi/(Tp*sqrt(1-zeta^2));

dynamics = table(wn, zeta)