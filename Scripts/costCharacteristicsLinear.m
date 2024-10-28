clc; clear; close all;

% Common for all plots
lWidth = 1.25;
costHeight = 21;

%% Motor(s)

% Cost [-]
wM = 2;
DMmax = 1000; % [cm^3/rev]
nMotor = 5;

motorSize = 1 : DMmax; % [cm^3/rev]
for idx = 1 : length(motorSize)
    costMotors(idx) = wM*(1 + motorSize(idx)/DMmax);
end

figure; hold on;
for i = 1 : nMotor
    col = [(1 - (i - 1)/nMotor) .5 .5];
    plot(motorSize,costMotors*i,'color',col, 'LineWidth', lWidth)
end
colororder('gem')
legend('1 Motor','2 Motors','3 Motors','4 Motors','5 Motors', ...
       'location','northwest', 'interpreter','latex')
title('Cost of \textbf{Motor(s)}', 'interpreter','latex')
xlabel('Motor Size $\left[\frac{cm^3}{rev}\right]$', 'interpreter','latex')
ylabel('Cost [-]', 'interpreter','latex')
ylim([0 costHeight])

%% Servo Valve(s)

% Cost [-]
wSV = 6;
QratMax = 1500; % [L/min]
nServo = 5;

servoSize = 1 : QratMax; % [L/min]
for idx = 1 : length(servoSize)
    costServos(idx) = wM*(1 + servoSize(idx)/QratMax);
end

figure; hold on;
for i = 1 : nServo
    col = [.5 .5 (1 - (i - 1)/nServo)];
    plot(servoSize,costServos*i,'color',col, 'LineWidth',lWidth)
end
legend('1 Servo','2 Servos','3 Servos','4 Servos','5 Servos', ...
       'location','northwest', 'interpreter','latex')
title('Cost of \textbf{Servo Valve(s)}', 'interpreter','latex')
xlabel('Servo Size $\left[\frac{L}{min}\right]$', 'interpreter','latex')
ylabel('Cost [-]', 'interpreter','latex')
ylim([0 costHeight])

%% Proportional Valve(s)

% Cost [-]
wPV = 4;
QnomMax = 1150; % [L/min]
nPvalve = 5;

pvalveSize = 1 : QnomMax; % [L/min]
for idx = 1 : length(pvalveSize)
    costPvalves(idx) = wM*(1 + pvalveSize(idx)/QratMax);
end

figure; hold on;
for i = 1 : 5
    col = [.5 (1 - (i - 1)/nPvalve) .5];
    plot(pvalveSize,costPvalves*i,'color',col, 'LineWidth',lWidth)
end
legend('1 Servo','2 Servos','3 Servos','4 Servos','5 Servos', ...
       'location','northwest', 'interpreter','latex')
title('Cost of \textbf{Proportional Valve(s)}', 'interpreter','latex')
xlabel('Servo Size $\left[\frac{L}{min}\right]$', 'interpreter','latex')
ylabel('Cost [-]', 'interpreter','latex')
ylim([0 costHeight])

%% Counterbalance Valve(s)

% Cost [-]
wCBV = 1;
QcapMax = 480; % [L/min]
nCBV = 5;

cbvType = 1 : QcapMax; % [L/min]
for idx = 1 : length(cbvType)
    costCBV(idx) = wM*(1 + cbvType(idx)/QratMax);
end

figure; hold on;
for i = 1 : nCBV
    col = [(1 - (i - 1)/nCBV) 0.25 (1 - (i - 1)/nCBV)];
    plot(costCBV*i,'color',col, 'LineWidth',lWidth)
end
legend('1 CBV','2 CBV','3 CBV','4 CBV','5 CBV', ...
       'location','northwest', 'interpreter','latex')
title('Cost of \textbf{Counterbalance Valve(s)}', 'interpreter','latex')
xlabel('Counterbalance Size $\left[\frac{L}{min}\right]$', 'interpreter','latex')
ylabel('Cost [-]', 'interpreter','latex')
ylim([0 costHeight])