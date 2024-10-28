clc; clear; close all;

%% Simulation as proof-of-concept: zRef subtraction

Tw = 10; % [sec]
zw_A = 1.2; % Amplitude
% 
% zL0 = 1; % [m]
% zp0 = 10; % [m]
% zSP = 5; % [m]

m = 24e3; % [kg]
g = 9.81; % [m/s^2]
dD = 0.45; % [m]
dR = 0.50; % [m]
dP = 0.15; % [m]
nSH = 3; % Number of Sheaves
nM = 5;  % Number of Motors
ig = 7;  % Gear Ratio
% 
% % From datasheet - Size 160
% J_motor = 0.022; % [kg*m^2]
% DM = 160.4e6; % [m^3/rev]


% star = (dD*dP) / (4*nSH*dR*ig*nM); % Total Transmission
% J = (1/2) * m * star; % [kg*m^2]

% zL = zL0; % [m]
% zLDot = 0; % [m/s]
% zLDotDot = 0;
% zw = zw_A*sin(0); % [m]

% % Controller gain
% Kp = 1; 

%%%%% Simulation %%%%%
time = 0;
step = 1e-5;
stopTime = 30;
idx = 1;
while time <= stopTime

    % if idx > 1
    %    zL = zL - zw; 
    % end
    % % Real example with amplitude and frequency
    zw = zw_A * sin(2*pi/Tw * time); % [m]
    % zwDot = 2*pi/Tw * zw_A * cos(2*pi/Tw * time); % [m/s]
    % zwDotDot = - 4*pi^2/(Tw^2) * zw_A * sin(2*pi/Tw * time); % [m/s^2]
    % zL = zL + zw;

    % Ideal example
    % zw = sin(time);
    % zwDot = cos(time);
    % zwDotDot = -sin(time)

    % zRef = zSP - zw; % [m]

    % % P-controller
    % e = zRef - zL;
    % MM = Kp*e;

    % zL = zL0 + zw; % [m]
    % % Physical Limits
    % if zL < 0
    %     zL = 0;
    % elseif zL > zp0
    %     zL = zp0;
    % end

    % zL_ext = zL + zw; % [m]

    plotTime(idx) = time; % [sec]
    % plotZref(idx) = zRef; % [m]
    plotZw(idx) = zw; % [m]
    % plotZwDot(idx) = zwDot; % [m/s]
    % plotZwDotDot(idx) = zwDotDot; % [m/s^2]
    % plotZL_ext(idx) = zL_ext;
    % plotZL(idx) = zL;

    % zLDotDot = (pL*DM/(2*pi) - MM) / (star*J); % [m/s^2]
    % Update
    % zLDot = zLDot + zLDotDot*step;
    % zL    = zL    + zLDot   *step;

    time = time + step;
    idx = idx + 1;
end
%%%%% Simulation %%%%%

% figure
% plot(plotTime, plotZw, 'b')
% hold on
% plot(plotTime, plotZwDot)
% plot(plotTime, plotZwDotDot)
% plot(plotTime, plotZref, 'k')
% xlabel('Time [seconds]', 'Interpreter','latex')
% lgd = legend('$z_w$', '$\dot{z}_w$','$\ddot{z}_w$', '$z_{ref}$', ...
%              'interpreter', 'latex', 'location', 'eastoutside');
% fontsize(lgd, 15, 'points')
% 
% figure
% plot(plotTime, plotFFuRef)
% hold on
% plot(plotTime, plotZref, 'k')
% title('Feedforward with unity gain $G=1$', 'interpreter', 'latex')
% xlabel('Time [seconds]', 'Interpreter','latex')
% lgd = legend('$u_{ref}$', '$z_{ref}$', 'interpreter', 'latex',...
%     'Location', 'eastoutside');
% fontsize(lgd, 15, 'points')

figure
plot(plotTime, plotZw, 'k')
% figure
% plot(plotTime, plotZL, 'b')
% hold on
% plot(plotTime, plotZref, 'k')
% yline(zp0, '--k')