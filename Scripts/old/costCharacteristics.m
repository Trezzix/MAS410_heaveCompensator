clc; clear; close all;

%% Motor(s)

% Cost [-]
wM = 2;
DMmax = 1000e-3; % [m^3/rev]
nMotor = 5;

motorType = [4.93 10.3 12 16 22.9 28.1 32 45.6 56.1 63 80.4 ...
             90 106.7 125 160.4 180 200 250 355 500 710 1000] * ...
             1e-3; % [cm^3/rev]
for idx = 1 : length(motorType)
    costMotors(idx) = wM*(1 + motorType(idx)/DMmax);
end

figure; hold on;
for i = 1 : nMotor
    col = [(1 - (i - 1)/nMotor) 0 0];
    plot(motorType, costMotors*i,'color',col, 'LineWidth', 1)
end
colororder('gem')
legend('1 Motor','2 Motors','3 Motors','4 Motors','5 Motors', ...
       'location','northwest')
title('Cost of Motor(s)')
xlabel('Motor Size #')
ylabel('Cost [-]')

%% Servo Valve(s)

% Cost [-]
wSV = 6;
QratMax = 1500 / 1000 / 60; % [m^3/sec]
nServo = 5;

servoType = [30 60 80 150 250 350 550 1000 1500] / 1000 / 60; % [m^3/sec]
for idx = 1 : length(servoType)
    costServos(idx) = wM*(1 + servoType(idx)/QratMax);
end

figure; hold on;
for i = 1 : nServo
    col = [0 0 (1 - (i - 1)/nMotor)];
    plot(servoType, costServos*i,'color',col, 'LineWidth',1)
end
legend('1 Servo','2 Servos','3 Servos','4 Servos','5 Servos', ...
       'location','northwest')
title('Cost of Servo Valve(s)')
xlabel('Servo Size #')
ylabel('Cost [-]')