clc;clear;close all;
%format bank
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
eta_hmM = 1; %hydromechanical effiency of motor
eta_vM = 0.94; %volumetric effiency of motor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chosen constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps = 220e5; %210 bar -> Pa
nm = 3; %number of motors
    %motor size 63
%Jm = 0.0012; %moment of inertia of motor
%Dm = 22.9e-6; %cm^3 -> m^3
    %proportional valve
npv = 1; %number of proportional valves
pr = 10e5; %bar -> Pa
%safety_factor = 1.1; % from lecture
deltaP_spool = 4.2e5; %bar -> pa, from datasheet of cvg50 31-08
deltaP_comp = 6e5; % bar -> pa, pg 12, cvg50
    %counter balance valve
ncbv = 2; %number of proportional valves
pM_in_lower = 18e5; %choose between 10-30 [bar], higher = bigger pressure in p1, lower = cavitation
pcr2_over = 1.3; %1.1 to 1.3
%Qr = 350 / (6*10^4);%L/min -> m^3/s
    %TODO: liquid
rou = 875; %placeholder
beta = 9e6; %very placeholder


%gear ratios
ip = (dR/2)/(dp/2);
iT = ig * ip; %total gear ratio between motor and drum
star = (dD/2)/(2*n_sh*ig*ip*nm); %[m]

%calculating max speed
zdot_max = (Zw * 2*pi) / Tw;
thetadot_D_max = (4*n_sh*zdot_max)/dD %max speed of drum [rad/s], absolute value, wrong?
thetadot_m_max = iT*thetadot_D_max %[rad/s] 
thetadot_m_max_rpm = thetadot_m_max * (60/(2*pi)) %rad/s -> rpm
M_M_max = ((mpl * g * dD * dp) / (4 * n_sh * dR * ig * nm)) * (1 + mu_eq * tanh(thetadot_m_max/w0))
%max acceleration
Zdotdot_max = -((Zw*(2*pi)^2)/(Tw^2)); %m/s^2
thetadotdot_D_max = -(4*n_sh*Zdotdot_max)/dD; %rad/s^2
thetadotdot_m_max = iT*thetadotdot_D_max; %[rad/s^2]

%for choosing motor size
pL_assume = (2/3) * ps; %chosen, changes for circuit B
D_min = (2*pi * M_M_max) / pL_assume; %ps -> delta_m_p(?) %TODO: might be missing volumetric efficiency
D_min_cm = D_min * 1e6
D_min_rpm = 63; %largest motor we can use with 
%auto choosing motor size
motorType = [4.93 10.3 12 16 22.9 28.1 32 45.6 56.1 63 80.4 ...
             90 106.7 125 160.4 180 200 250 355 500 710 1000];
MotorJ = [0.00006 0.0004 0.0004 0.0004 0.0012 0.0012 0.0012 0.0024 0.0042 0.0042 0.0072 0.0072 ...
    0.0116 0.0116 0.0220 0.0220]; %dosent include all of them
if D_min_cm > D_min_rpm
    error("Motor size greater than 63, increase nm")
end

for i_for = 1:length(motorType)
    if motorType(i_for) > D_min_cm
        Dm_cm = motorType(i_for)
        Dm = motorType(i_for) * 1e-6;
        Jm = MotorJ(i_for);
        break
    end
end
Jtot = Jm + (1/2)*(mpl)*star^2; %double check

%pressure
pL_max = (M_M_max + Jtot * thetadotdot_m_max) * ((2*pi)/Dm);
%pL = (M_M_max) * ((2*pi)/Dm);
pL_max_bar = pL_max*1e-5
%flow
Qm_max = eta_vM * (Dm/(2*pi)) * thetadot_m_max;
Qm_max_Lpmin = Qm_max * 6*10^4 %for prop valve sizing

%%%proportional valve - spool
%type: closed center, symmetric, CVGxx 31-xx
Qm_max_total = (Qm_max*nm)/(npv);
Qm_max_total_lMin = Qm_max_total * 6*10^4;
spoolFlows = [60, 140, 220, 500, 820, 1000, 950, 1150];
spoolTypes = ["CVG30 31-00", "CVG30 31-01", "CVG30 31-02", "CVG30 31-05",...
    "CVG50 31-08", "CVG50 31-10", "CVG60 31-10", "CVG60 31-20"];

for i_for = 1:length(spoolFlows)
    if spoolFlows(i_for) > Qm_max_total_lMin
        Q_nom_spool = spoolFlows(i_for)
        chosenSpool = spoolTypes(i_for)
        break
    end
end
%crank
%Qr_spool = 1550 / 6e4; %from datasheet of CVG50 31-10
Qr_spool = 1220 / 6e4; %from datasheet of CVG50 31-08
Qr_comp = 1150 / 6e4; %from datasheet of CVG50 pg 12

CdAd_spool = Qr_spool/sqrt((2/rou) * pr); %m^3 main spool
%CdAd_comp = Qr_comp/sqrt((2/rou) * pr); %m^3 can this be done? - doubt since Ad changes over time

pcr1 = (Qm_max_total^2 * rou) / (CdAd_spool^2 *2); %3.36

%pcr1 = rou/2 * (Qm_max_total)^2 * (1/(CdAd_spool^2) + 1/(CdAd_comp^2)) - 6e5;

%delta_p_spool = (Qm_max_total^2 * rou) / (CdAd_spool^2 * 2);
pcr1_bar = pcr1 * 1e-5
screw_flow_open = 450;
screw_flow_max = 1300;
screw_max_turns = 9.5;

pcr1_turns = pcr1_bar/(11/screw_max_turns) %TODO: double check, 11??? hvor kommer det fra..

pM_in = ps - (deltaP_comp + deltaP_spool);

%%%counterbalance / overbalance valve:
cbv_alpha_list = [1.5 2 2.3 3 4.5 10];
cbv_name_list = ["CBiB" "CBIY" "CBIL" "CBIA" "CBIG" "CBIH"]
pM_out = pM_in - pL_max; %unsure if correct, TODO: doublec check
pcr2 = pL_max * pcr2_over; % pL_max is very close to the method used in the examples
pcr2_bar = pcr2 * 1e-5
pRet = (Qm_max_total^2 * rou)/(CdAd_spool^2 * 2);
pRet_bar = pRet * 1e-5
p1 = ((M_M_max * 2 * pi) / Dm);
alpha_max = (p1 +pM_in_lower - pcr2 - pRet) / (pRet - pM_in_lower); % wrong?

for i_for = 1:length(cbv_alpha_list)
    if cbv_alpha_list(i_for) > alpha_max
        alpha_cbv = cbv_alpha_list(i_for-1)
        cbv_type = cbv_name_list(i_for-1)
        break
    end
end

pM_in_lower = (p1 - pcr2 + pRet*(-1 -alpha_cbv))/((-1 -alpha_cbv));
pM_in_lower_bar = pM_in_lower * 1e-5
max_capacity = 480; %L/min
n_cbv_min = ceil(Qm_max_total_lMin/max_capacity)