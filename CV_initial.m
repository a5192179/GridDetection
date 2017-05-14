function [X,Xpda,P_initial] = CV_initial(target_meas)
%   This function is applied to initialize the ture track
%   The Model Probabilities are adjustable here
global R T total_frame
X = zeros(4,total_frame);
X(:,1) = [target_meas(1,1);(target_meas(1,2)-target_meas(1,1))/T;...
          target_meas(2,1);(target_meas(2,2)-target_meas(2,1))/T];
X(:,2) = [target_meas(1,2);(target_meas(1,2)-target_meas(1,1))/T;...
          target_meas(2,2);(target_meas(2,2)-target_meas(2,1))/T];
P_initial   = [R(1,1),R(1,1)/T,0,0;
               R(1,1)/T,2*R(1,1)/T^2,0,0;
               0,0,R(2,2),R(2,2)/T;
               0,0,R(2,2)/T,2*R(2,2)/T^2];
Xpda = X(:,2);