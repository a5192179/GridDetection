clc
clear all
close all
ss = 4; % state size
os = 2; % observation size
F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
H = [1 0 0 0; 0 1 0 0];
Q = 0.1*eye(ss);
R = 50*eye(os);
initx = [52000 2000 -50 -50]';
initV = 10*eye(ss);
x=[(52000:-50:48000)',(2000:-50:-2000)']';
y=x+100*rand(size(x));
[xfilt, Vfilt, VVfilt, loglik] = kalman_filter(y, F, H, Q, R, initx, initV);
figure(1)
plot(x(1,:), x(2,:), 'ks-');
hold on
plot(y(1,:), y(2,:), 'g*');
hold on
plot(xfilt(1,:), xfilt(2,:), 'rx:');
legend('true', 'observed', 'filtered', 3)