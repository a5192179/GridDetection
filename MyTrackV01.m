function MyTrack(Area,Target_Location,Velocity,gridWidth,ResultSaveFolder,mont_num)
global total_frame T Pfa Pd Pg SNR lamda gama

%% Basic Parameters
% saveRootFolder='..\..\Data';
total_frame = length(Target_Location(:,1));
T= 1; % the interval of scanning
% total_time = total_frame*T; % the total frame of scanning
SNR = 15; % Sigal to Noise Ratio
Pfa = 0.005; % the expected false alarm rate
lamda = Pfa; % clutter density that has relationship with Pfa (For normalizing the resolution cells, it's equal to Pfa)
Pg = 0.9997; % the probability of the measurement that fall in the validation region
gama = 4; % gate size parameter
% R = diag([1,1]);
% mont_num = 10; % Monte carlo number
% Save_data_detection=zeros(250,250,5);%���޸ġ�20150124-��˧-����ÿ�β�����data_detection
% Target_Location=[190,0;180,0;170,0;160,0;150,0];
% Velocity=[-10,0];
%% CFAR Parameters
% N_Azi = 7; N_APro = 1; N_Range = 7; N_RPro = 1; % the length of sliding window and protect cells both at azimuth and range
% N_aver = N_Azi*N_Range - (2*N_APro + 1)*(2*N_RPro + 1); % the number of auxiliary cells
% alpha = N_aver*(Pfa^(-1/N_aver) - 1); % threshold product factor ��Fundamental of Radar Sigal Processing��p.265
Pd = 1e-5;%(1+alpha/(N_aver*(1+10^(SNR/10))))^(-N_aver); % detection probability ��Fundamental of Radar Sigal Processing��p.266

%% Rayleigh distribution Parameters
% sigma = 1; % Rayleigh distribution, average power (Normalization)

%% Calculate RMSE, some matrices
T_rmse =zeros(total_frame,mont_num);
T_mse = zeros(1,total_frame);
Loss = 0; % Track Loss flag
C1 = gridWidth^2; C2 = 2*gridWidth^2; % Track Loss Threashood (lower,upper)
success = 0; % flag of tracks which are successfully tracked
mode = 0; % if mode == 1,using AI,if mode == 0, common mode

%% Monte Carlo Test
for i_Montecarlo = 1:1%mont_num
    %%
    [target_state,target_meas] = non_maneuvering; % the target location������
    target_state=[Target_Location(:,1)';repmat(Velocity(1,1),1,total_frame);Target_Location(:,2)';repmat(Velocity(1,2),1,total_frame)];
    target_meas=[Target_Location(:,1)';Target_Location(:,2)'];
%     [target_state,target_meas] = non_maneuvering; % the target location������
    [X,Xpda,P_estimate] = CV_initial(target_meas); % Initialization--PDAF parameters (the first two frames)������
%     data_target = raylrnd(sqrt(sigma+10^(SNR/10)),1,total_frame); % the target return (Rayleigh)
    k = 0; l = 0;
    
    %% Estimate the target state, output: trajectory
    for jj = 3:total_frame
        %% Generate random data
        %           data = raylrnd(sqrt(sigma),x_num,y_num); % Rayleigh random noise data per frame
        %           data(target_meas(1,jj),target_meas(2,jj)) = data_target(jj); % Data:target pluse noise
        
        %% Track maintance
        %           [data_detection] = CA_CFAR(data,N_Azi,N_APro,N_Range,N_RPro,alpha); % CFAR processing
        %           Save_data_detection(:,:,jj)=data_detection;
        detectionResultSaveFileName=[ResultSaveFolder,'\',num2str(jj,'%02d'),'.txt'];
        Result_Detection=importdata(detectionResultSaveFileName);
        [Xpda,P_estimate] = PDAF(Result_Detection,Xpda,P_estimate,mode); % PDAF������
        %           [Xpda,P_estimate] = PDAF(data_detection,Xpda,P_estimate,mode); % PDAF������
        X(:,jj) = Xpda; % target state estimation (storage)
        
        
        %% Trajectory quality evaluation
        T_mse(:,jj) = (X(1,jj)-target_state(1,jj))^2 + (X(3,jj)-target_state(3,jj))^2; % square error per frame
        if sqrt(T_mse(:,jj)) > C1 & sqrt(T_mse(:,jj)) < C2
            k = k+1;
        elseif sqrt(T_mse(:,jj))<C1
            k = 0;
        elseif sqrt(T_mse(:,jj)) > C2
            l = l + 1;
        end
        if k == 8
            k = 0; l = l + 1;
        end
    end
    if l == 0
        T_rmse(:,i_Montecarlo) = T_mse;%ÿ�з�һ�����ؿ����֡�����
        % MIU(:,:,ii) = Miu; %���ף��ң�
        success = success + 1;
    else
        Loss = Loss + 1;
    end
    
end

%% Statisticsͳ��
%[clutter_x,clutter_y] = find(data_detection); % Radar clutter plots
%���޸�1����˧-20150124-ԭ���������ƽ����ѭ���⣨�����ã��Ҹĵ���ͼ������
RMSE = zeros(1,total_frame);
% % Miu_aver = zeros(3,total_frame); %���ף�
%Miu_aver = zeros(2,total_frame); %���ң�
for i_frame = 1:total_frame
    RMSE(:,i_frame) = sqrt(sum(T_rmse(i_frame,:))/success);
    % Miu_aver(1,ii) = sum(MIU(1,ii,:))/success; % ���ף��ң�
    %Miu_aver(2,ii) = sum(MIU(2,ii,:))/success; % ���ף��ң�
    % %     Miu_aver(3,ii) = sum(MIU(3,ii,:))/success; % ���ף�
end

%% Display the results
success
Loss
% figure(1)
% plot(clutter_x,clutter_y,'r*')
% hold on
% plot(X(1,:),X(3,:),'k-',target_state(1,:),target_state(3,:))
% title(['Trajectory Display  SNR=',num2str(SNR),'']);
% xlabel('X-coordinate');
% ylabel('Y-coordinate');
% legend('Clutter','Estimate','True');
%
% figure(2)
% plot(RMSE,'k')
% title('The Root Mean Square Error Curve');
% xlabel('Frame');
% ylabel('RMSE');
figure
for i=1:total_frame
    if i>2
        %         [clutter_x,clutter_y] = find(Save_data_detection(:,:,i));% Radar clutter plots
        detectionResultSaveFileName=[ResultSaveFolder,'\',num2str(i,'%02d'),'.txt'];
        Result_Detection=importdata(detectionResultSaveFileName);
        plot(Result_Detection(:,1),Result_Detection(:,2),'m*')
    end
    hold on
    if i==1
        plot(X(1,i),X(3,i),'b.-',target_state(1,i),target_state(3,i),'ro-')
    else
        plot(X(1,1:i),X(3,1:i),'b.-',target_state(1,1:i),target_state(3,1:i),'ro-')
    end
    hold off
    title(['Trajectory Display  SNR=',num2str(SNR),'']);
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    grid on
    axis([Area(1,1) Area(1,2) Area(2,1) Area(2,2)])
    legend('Detect','Estimate','True');
    print(gcf,'-dbmp',strcat('../Figure/lsmycv',num2str(i),'.bmp'));
end
filename1='../Figure/lsmycv.gif';

for iter=1:total_frame
    im=imread(['../Figure/lsmycv',num2str(iter),'.bmp']);
    [imind,map2] = rgb2ind(im,256);
    if iter==1
        imwrite(imind,map2,filename1,'gif', 'Loopcount',inf,'DelayTime',0.2)
    else
        imwrite(imind,map2,filename1,'gif','writeMode','append','delaytime',0.2);
    end
end

figure
plot(RMSE,'k')
title('The Root Mean Square Error Curve');
xlabel('Frame');
ylabel('RMSE');


