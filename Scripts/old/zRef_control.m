clc; clear; close all;

%% Simulation as proof-of-concept: zRef subtraction

Tw = 10; % [sec]
zw_A = 1.2; % Amplitude

zL0 = 1; % [m]
zp = 10; % [m]
% m = 24e3; % [kg]
m = 1;
g = 9.81; % [m/s^2]

zL = zL0; % [m]
zLDot = 0; % [m/s]
zw = zw_A*sin(0); % [m]

% Controller gain
Kp = 1; 

%%%%% Simulation %%%%%
time = 0;
step = 1e-5;
stopTime = 30;
idx = 1;
while time <= stopTime

    if idx > 1
       zL = zL - zw; 
    end
    % Real example with amplitude and frequency
    zw = zw_A * sin(2*pi/Tw * time); % [m]
    zwDot = 2*pi/Tw * zw_A * cos(2*pi/Tw * time); % [m/s]
    zwDotDot = - 4*pi^2/(Tw^2) * zw_A * sin(2*pi/Tw * time); % [m/s^2]
    zL = zL + zw;

    % Ideal example
    % zw = sin(time);
    % zwDot = cos(time);
    % zwDotDot = -sin(time)

    zRef = zp - zw; % [m]

    % P-controller
    e = zRef - zL;
    F = Kp*e;

    % zL = zL0 + zw; % [m]
    % Physical Limits
    if zL < 0
        zL = 0;
    elseif zL > zp
        zL = zp;
    end

    % Feedforward testing
    zRefDot = - zwDot;
    FFuRef = zRef + zRefDot; % Assuming unity gain

    plotTime(idx) = time; % [sec]
    plotZref(idx) = zRef; % [m]
    plotZw(idx) = zw; % [m]
    plotZwDot(idx) = zwDot; % [m/s]
    plotZwDotDot(idx) = zwDotDot; % [m/s^2]
    plotFFuRef(idx) = FFuRef;
    plotZL(idx) = zL;

    zLDotDot = (F + m*g) / m; % [m/s^2]
    % Update
    zLDot = zLDot + zLDotDot*step;
    zL    = zL    + zLDot   *step;

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
plot(plotTime, plotZL, 'b')
hold on
plot(plotTime, plotZref, 'k')
yline(zp, '--k')