 function [target_state,target_meas] = non_maneuvering
% This funcion is applied to generate the target's track
% target_state: the target's state vecter, which has 4 varibles-the x,y
%               positions and x,y velocities
% target_meas : the radar measurements of the target, it includes two
%               elements-x,y location, and both of the elements must be
%               integer
% close all
% clear all
% clc

global total_frame T F1 H Q1 R
% total_time = 50;
T = 1;
% total_frame = total_time/T;
F1 = [1,T,0,0;
      0,1,0,0;
      0,0,1,T;
      0,0,0,1]; 
H = [1,0,0,0;
     0,0,1,0];
q = 0.001; % process noise power spectrum parameter
Q1 = q*[T^3/3,T^2/2,0,0;
        T^2/2,T,0,0;
        0,0,T^3/3,T^2/2;
        0,0,T^2/2,T];
R = diag([1,1]);

target_state(:,1)=[5;4;5;4]; % initial position
measure_noise = sample_gaussian(zeros(length(R),1),R,total_frame)';
target_meas(:,1)= H*target_state(:,1)+measure_noise(:,1);
process_noise = sample_gaussian(zeros(length(Q1),1),Q1,total_frame)';

for ii=2:total_frame
     target_state(:,ii) = F1*target_state(:,ii-1) + process_noise(:,ii-1);
     target_meas(:,ii) = H*target_state(:,ii) + measure_noise(:,ii);
end

target_meas = round(target_meas);
 
 
% plot(target_state(1,:),target_state(3,:),'k')
% hold on
% plot(target_meas(1,:),target_meas(2,:),'r*')
% % axis([0 150 0 150]);
% grid
