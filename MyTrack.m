function [l,T_mse,Xpda]=MyTrack(Area,Target_Location,Velocity,gridWidth,Result_Detection,flag_FigureTrack)
global total_frame T Pfa Pd Pg SNR lamda gama i_Figure
%% Basic Parameters
% saveRootFolder='..\..\Data';
total_frame = length(Target_Location(:,1));
T= 1; % the interval of scanning
% total_time = total_frame*T; % the total frame of scanning
SNR = 10; % Sigal to Noise Ratio
Pfa = 1e-5;%0.005; % the expected false alarm rate
lamda = Pfa; % clutter density that has relationship with Pfa (For normalizing the resolution cells, it's equal to Pfa)
Pg = 0.9997; % the probability of the measurement that fall in the validation region
gama = 4; % gate size parameter
% R = diag([1,1]);
% mont_num = 10; % Monte carlo number
% Save_data_detection=zeros(250,250,5);%【修改】20150124-李帅-保存每次产生的data_detection
% Target_Location=[190,0;180,0;170,0;160,0;150,0];
% Velocity=[-10,0];
%% CFAR Parameters
% N_Azi = 7; N_APro = 1; N_Range = 7; N_RPro = 1; % the length of sliding window and protect cells both at azimuth and range
% N_aver = N_Azi*N_Range - (2*N_APro + 1)*(2*N_RPro + 1); % the number of auxiliary cells
% alpha = N_aver*(Pfa^(-1/N_aver) - 1); % threshold product factor 《Fundamental of Radar Sigal Processing》p.265
Pd = 1;%(1+alpha/(N_aver*(1+10^(SNR/10))))^(-N_aver); % detection probability 《Fundamental of Radar Sigal Processing》p.266

%% Rayleigh distribution Parameters
% sigma = 1; % Rayleigh distribution, average power (Normalization)

%% Calculate RMSE, some matrices
% T_rmse =zeros(total_frame,mont_num);
T_mse = zeros(1,total_frame);
% Loss = 0; % Track Loss flag
C1 = gridWidth^2; C2 = 2*gridWidth^2;%C1 = 5; C2 = 20; % Track Loss Threashood (lower,upper)
% success = 0; % flag of tracks which are successfully tracked
mode = 0; % if mode == 1,using AI,if mode == 0, common mode

%% Monte Carlo Test
% for i_Montecarlo = 1:1%mont_num
    %%
    [target_state,target_meas] = non_maneuvering(gridWidth); % the target location（丙）
    target_state=[Target_Location(:,1)';repmat(Velocity(1,1),1,total_frame);Target_Location(:,2)';repmat(Velocity(1,2),1,total_frame)];
    target_meas=[Target_Location(:,1)';Target_Location(:,2)'];
%     [target_state,target_meas] = non_maneuvering; % the target location（丙）
    [X,Xpda,P_estimate] = CV_initial(target_meas); % Initialization--PDAF parameters (the first two frames)（丙）
%     data_target = raylrnd(sqrt(sigma+10^(SNR/10)),1,total_frame); % the target return (Rayleigh)
    k = 0; l = 0;
    
    %% Estimate the target state, output: trajectory
    num_Inside=0;
    for i_Frame = 3:total_frame
        %% Generate random data
        %           data = raylrnd(sqrt(sigma),x_num,y_num); % Rayleigh random noise data per frame
        %           data(target_meas(1,jj),target_meas(2,jj)) = data_target(jj); % Data:target pluse noise
        
        %% Track maintance
        %           [data_detection] = CA_CFAR(data,N_Azi,N_APro,N_Range,N_RPro,alpha); % CFAR processing
        %           Save_data_detection(:,:,jj)=data_detection;
%         detectionResultSaveFileName=[ResultSaveFolder,'\',num2str(i_Frame,'%02d'),'.txt'];
%         Result_Detection=importdata(detectionResultSaveFileName);
        Result_Detection_ThisFrame=zeros(100,2);
        Result_Detection_ThisFrame(1:100,1:2)=Result_Detection(i_Frame,:,:);
        Result_Detection_ThisFrame(Result_Detection_ThisFrame==1.211)=[];
        if isempty(Result_Detection_ThisFrame)
            
        end
        [Xpda,P_estimate,flag_Inside] = PDAF(Result_Detection_ThisFrame,Xpda,P_estimate,mode); % PDAF（丙）
        num_Inside=num_Inside+flag_Inside;
        %           [Xpda,P_estimate] = PDAF(data_detection,Xpda,P_estimate,mode); % PDAF（丙）
        X(:,i_Frame) = Xpda; % target state estimation (storage)
        
        
        %% Trajectory quality evaluation
        T_mse(:,i_Frame) = (X(1,i_Frame)-target_state(1,i_Frame))^2 + (X(3,i_Frame)-target_state(3,i_Frame))^2; % square error per frame
        if sqrt(T_mse(:,i_Frame)) > C1 & sqrt(T_mse(:,i_Frame)) < C2
            k = k+1;
        elseif sqrt(T_mse(:,i_Frame))<C1
            k = 0;
        elseif sqrt(T_mse(:,i_Frame)) > C2
            l = l + 1;
        end
        if k == 8
            k = 0; l = l + 1;
        end
    end
%     if l == 0
%         T_rmse(:,i_Montecarlo) = T_mse;%每列放一次蒙特卡洛各帧的误差
%         % MIU(:,:,ii) = Miu; %（甲，乙）
%         success = success + 1;
%     else
%         Loss = Loss + 1;
%     end
    
% end



%% Display the results
% success
% Loss
% figure()
% plot(clutter_x,clutter_y,'r*')
% hold on
% plot(X(1,:),X(3,:),'k-',target_state(1,:),target_state(3,:))
% title(['Trajectory Display  SNR=',num2str(SNR),'']);
% xlabel('X-coordinate');
% ylabel('Y-coordinate');
% legend('Clutter','Estimate','True');

% figure(2)
% plot(RMSE,'k')
% title('The Root Mean Square Error Curve');
% xlabel('Frame');
% ylabel('RMSE');
%%
if flag_FigureTrack==1
    figure_TrackPoint=figure(i_Figure);
    i_Figure=i_Figure+1;
    figure(figure_TrackPoint)
    for i_Frame=1:total_frame
        plot_True=plot(target_state(1,1:i_Frame),target_state(3,1:i_Frame),'ro-');
        %     legend(plot_True,'True')
        hold on
        %     if i>2
        %         [clutter_x,clutter_y] = find(Save_data_detection(:,:,i));% Radar clutter plots
%         detectionResultSaveFileName=[ResultSaveFolder,'\',num2str(i_Frame,'%02d'),'.txt'];
%         Result_Detection=importdata(detectionResultSaveFileName);
        Result_Detection_ThisFrame=Result_Detection(i_Frame,:,:);
        Result_Detection_ThisFrame(Result_Detection_ThisFrame==1.211)=[];
        plot_Detect=plot(Result_Detection_ThisFrame(:,1),Result_Detection_ThisFrame(:,2),'m*');
        %     legend(plot_Detect,'Detect')
        hold on
        %     end
        
        plot_Estimate=plot(X(1,1:i_Frame),X(3,1:i_Frame),'b.-');
        legend([plot_True,plot_Detect,plot_Estimate],'True','Detect','Estimate','Location','NorthWest')
        hold on
        
        %     hold off
        title(['跟踪轨迹 SNR=',num2str(SNR),'dB 栅格宽度=',num2str(gridWidth),'m']);
        xlabel('X/m');
        ylabel('Y/m');
        grid on
        axis([Area(1,1) Area(1,2) Area(2,1) Area(2,2)])
        %     legend('Estimate','True');
        print(gcf,'-dbmp',strcat('../Figure/lsmycv',num2str(i_Frame),'.bmp'));
    end
end
% filename1='../Figure/lsmycv.gif';
% 
% for iter=1:total_frame
%     im=imread(['../Figure/lsmycv',num2str(iter),'.bmp']);
%     [imind,map2] = rgb2ind(im,256);
%     if iter==1
%         imwrite(imind,map2,filename1,'gif', 'Loopcount',inf,'DelayTime',0.2)
%     else
%         imwrite(imind,map2,filename1,'gif','writeMode','append','delaytime',0.2);
%     end
% end
%%
% figure
% plot(RMSE,'k')
% title('The Root Mean Square Error Curve');
% xlabel('Frame');
% ylabel('RMSE');


