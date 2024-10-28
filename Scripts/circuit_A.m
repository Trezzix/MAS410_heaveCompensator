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
nm = 1; %number of motors
    %motor size 63
%Jm = 0.0012; %moment of inertia of motor
%Dm = 22.9e-6; %cm^3 -> m^3
    %valve
nv = 1; %number of valves
pr = 10e5; %bar -> Pa ,datasheet
safety_factor = 1.1; % from lecture
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
pL_assume = (2/3) * ps; %chosen
D_min = (2*pi * M_M_max) / pL_assume; %ps -> delta_m_p(?) %TODO: might be missing volumetric efficiency
D_min_cm = D_min * 1e6
D_min_rpm = 80; %largest motor we can use with 
%auto choosing motor size
motorType = [4.93 10.3 12 16 22.9 28.1 32 45.6 56.1 63 80.4 ...
             90 106.7 125 160.4 180 200 250 355 500 710 1000];
MotorJ = [0.00006 0.0004 0.0004 0.0004 0.0012 0.0012 0.0012 0.0024 0.0042 0.0042 0.0072 0.0072 ...
    0.0116 0.0116 0.0220 0.0220 0.0353 0.061 0.102 0.178 0.55 0.55]; %kg/m^2
if D_min_cm > D_min_rpm
    warning("Motor size greater than 80, increase nm")
end

for i_for = 1:length(motorType)
    if motorType(i_for) > D_min_cm
        Dm_cm = motorType(i_for)
        Dm = motorType(i_for) * 1e-6;
        Jm = MotorJ(i_for);
        break
    end
end
Jtot = Jm + (mpl)*star^2; %TODO: double check

%pressure
pL_max = (M_M_max + Jtot * thetadotdot_m_max) * ((2*pi)/Dm);
%pL = (M_M_max) * ((2*pi)/Dm);
pL_max_bar = pL_max*1e-5
%flow
Qm_t = (Dm/(2*pi)) * thetadot_m_max;
Qm_t_Lpmin = Qm_t * 6*10^4
%no load method
%Qm_L_min = Qm_t * 6*10^4
Qm_NL = Qm_t * sqrt(ps/(ps - pL_max));
Qm_NL_Lpmin = Qm_NL * 6*10^4

%servo valve sizing
Qm_NL_total = (Qm_NL*nm)/(nv);
Qr_min = safety_factor * Qm_NL_total * sqrt(pr/ps);
Qr_min_lpmin = Qr_min * 6*10^4 % [l/min]
%auto choosing valve size
servoType = [30 60 80 150 250 350 550 1000 1500]; % [l/min]
for i_for = 1:length(servoType)
if servoType(i_for) > Qr_min_lpmin
    Qr_lpmin = servoType(i_for)
    Qr = servoType(i_for) / (6*10^4); %L/min -> m^3/s
    break
end
end

CdAd = Qr/sqrt((2/rou) * pr); %m^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = 0; %s
timeend = 10; %s
timestep = 1e-5; %s
counter = 1;

z= 0;
theta = 0;
thetadot = 0;
thetadotdto = 0;
zPl = 0;
p1 = 0; %TODO needs to be eq
p2 = 0;
p1dot = 0;
p2dot = 0;
Q1 = 0;
Q2 = 0;
V1 = 1e-6; %m^3, temporary, = 1 L
V2 = 1e-6; %m^3
V3 = 1e-6; %m^3
u = 1; % for now

tic
while time<timeend
    %simulate wave
    z = Zw * sin(((2*pi)/Tw)*time);
    zdot = Zw * cos((2*pi/Tw)*time) * ((2*pi)/Tw);
    zdotdot = ((-Zw*(2*pi)^2) / (Tw^2))* sin(((2*pi)/Tw) * time);
    %flow
    Q1 = CdAd*u*sign(ps - p1)*sqrt((2/rou) * abs(ps - p1));
    Q2 = CdAd*u*sign(p2)*sqrt((2/rou) * abs(p2));
    Qm = thetadot * Dm; %%%%%%%%%%%%%%?
    %diff eqs
    p1dot = (beta/V2) * (Q1 - Qm); %?
    p2dot = (beta/V3) * (Qm - Q2);
    PL_sim = p1 - p2;
    %M_M_sim = (mpl * g) / star;
    thetadotdot = (((PL_sim * Dm) / 2 * pi) - M_M_max)/Jtot; %TODO: M_M will change over time
    zPl = z + (theta * star); % this is correct

%save values
    z_graph(counter) = z;
    zpl_graph(counter) = zPl;
    PL_graph(counter) = PL_sim * 1e-5; % to bar
    P1_graph(counter) = p1 * 1e-5; % to bar
    P2_graph(counter) = p2 * 1e-5; % to bar
    Q1_graph(counter) = Q1;
    Q2_graph(counter) = Q2;
    time_graph(counter) = time;

    p1 = p1 + p1dot * timestep;
    p2 = p2 + p2dot * timestep;
    thetadot = thetadot + thetadotdot * timestep;
    theta = theta + thetadotdot * timestep;

    counter = counter + 1;
    time = time + timestep;

end
toc

figure
subplot(2,2,1)
plot(time_graph,z_graph)
hold on
plot(time_graph,zpl_graph)
grid on
legend('z-platform','z-payload')

subplot(2,2,2)
plot(time_graph,PL_graph)
hold on
plot(time_graph,P1_graph)
plot(time_graph,P2_graph)
grid on
legend('PL','P1','P2')

subplot(2,2,3)
plot(time_graph,Q1_graph)
hold on 
plot(time_graph,Q2_graph)
grid on
legend('Q1','Q2')