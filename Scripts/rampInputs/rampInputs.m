clc; clear; close all;

lWidth = 1.5;

% Initial values
t = -0.5; % [sec]
step = 1e-3; % [sec]
stopTime = 1.5; % [sec]
tR = 1;
idx = 1;

tic
while t <= stopTime    
    % Ramp Input: Linear
    if t < 0
        linear = 0;
    elseif t <= tR
        linear = (1 - 0) / (tR - 0) * t;
    else
        linear = 1;
    end
    
    % Ramp Input: Cubic polynomial
    tau = t/tR; % [-]
    if tau < 0
        smooth = 0;
    elseif tau <= 1
        smooth = 3*tau^2 - 2*tau^3;
    else
        smooth = 1;
    end

    % Save current iteration
    plotTime(idx) = t;
    plotSmooth(idx) = smooth;
    plotLinear(idx) = linear;
    
    % Increment
    t = t + step;
    idx = idx + 1;
end
toc

figure(Name='Linear')
plot(plotTime, plotLinear, 'k','LineWidth',lWidth)
title('Linear')
ylim([-0.1 1.1])
grid on

figure(Name='Cubic Polynomial')
plot(plotTime, plotSmooth, 'k','LineWidth',lWidth)
title('Cubic Polynomial')
ylim([-0.1 1.1])
grid on