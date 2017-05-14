%**************************************************************************
%*��Creat time����2017-01-27 15:08          ��Version����0.0
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*��������������������
%*
%*
%*
%*��Description����
%*
%*
%*-------------------------------------------------------------------------
%*��Modification����****-**-** **��**       ��Version����*.*
%*
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*
%*
%*
%*��Description����
%*
%*
%*
%**************************************************************************
clc
close all
clear all
fprintf('Begin time = %s\n',datestr(now))
%%
%��������
saveRootFolder='..\Data';
saveFolder_AlignmentResult=[saveRootFolder,'\DetectionResult\Alignment'];
if ~exist(saveFolder_AlignmentResult,'dir') 
    mkdir(saveFolder_AlignmentResult)
end 
saveFolder_ElectResult=[saveRootFolder,'\DetectionResult\Elect'];
if ~exist(saveFolder_ElectResult,'dir') 
    mkdir(saveFolder_ElectResult)
end 
method=2;%1���룬2��ѡ��3����
flag_Track=1;
flag_FigureDetect=1;
flag_FigureTrack=0;
%�״�λ�ñ�ţ�դ���С��ţ�������
order_Radar=4;
order_GridWidth=4;
order_Aera=1;
[Coord_T,GridWidth,Area]=SelectParameter(order_Radar,order_GridWidth,order_Aera);
%�״����
% Coord_T=[
%     0 50;%�����״�����1 2 3 4 5 XY
%     0 25;
%     0 0;
%     0 -25;
%     0 -50]*1e3;
Coord_R=Coord_T;
distanceUnit=50;%��λ��m
detectRange=300e3;%�״���Զ̽�����
num_distanceUnit=ceil(detectRange/distanceUnit);
%���泡��
% Area=[49.5,50.5;          %X
%     -0.5,0.5]*1e3;      %Y
%
Coord_TargetBegin=[50.4,0.4]*1e3;
% Coord_TargetBegin=[50.2,0.2]*1e3;
%
% Coord_target=Coord_targetBegin;
Velocity=[-10,-10];
Time=(0:79)';%79
num_Frame=length(Time);
% Time=0;
% TimeBegin=zeros(size(Time));
% Coord_target(:,1)=Coord_target(:,1)+velocity*(Time-TimeBegin)
% GirdWidth=[30;40;50;60;70;100;150];%��λ
% GridWidth=50;
%�������
% SNR=[linspace(-10,-2,5)';linspace(0,15,31)';linspace(16,20,5)'];
SNR=10;
num_Montecarlo=1;
%%
%��������
num_T=length(Coord_T(:,1));
num_Target=length(Coord_TargetBegin(:,1));
num_R=length(Coord_R(:,1));
num_Channel=num_T*num_R;
num_GridWidth=length(GridWidth(:,1));
global i_Figure
i_Figure=1;
%%
%��������
saveFileName_AlignmentGate=[saveRootFolder,'\Gate\AlignmentGate',num2str(num_Channel,'%02d'),'.txt'];
saveFileName_ElectGate=[saveRootFolder,'\Gate\ElectGate',num2str(num_Channel,'%02d'),'.txt'];
% if exist(saveFileName_AlignmentGate,'file')
%     gate_Alignment=importdata(saveFileName_AlignmentGate);
% else
%     gate_Alignment=importdata(indexSaveFileName);
% end
gate_Alignment=importdata(saveFileName_AlignmentGate);
% Gate_Elect=importdata(saveFileName_ElectGate);
%%
%��������
%%
%--------------------------------------------------------------------------
%���ٻ�ͼ����
if flag_Track==1
    figure_RMSE=figure(i_Figure);
    i_Figure=i_Figure+1;
    Mark_RMSE={'b-';'r:';'m-.';'k--';'g-'};
    Legend_RMSE=cell(num_GridWidth,1);
    RMSE = zeros(length(GridWidth(:,1)),num_Frame);
end
%--------------------------------------------------------------------------

for i_GridWidth=1:num_GridWidth
    gridWidth=GridWidth(i_GridWidth,1);%��λ
    %����դ��
    Temp_XCoord=(Area(1,1):gridWidth:Area(1,2))';%������
    num_XCoord=length(Temp_XCoord(:,1));
    Temp_YCoord=(Area(2,2):-gridWidth:Area(2,1))';%���ϵ���
    num_YCoord=length(Temp_YCoord(:,1));
    num_Grid=num_XCoord*num_YCoord;
    %----------------------------------------------------------------------
    %������������
    %�ж��ļ��Ƿ����
    if method==1||method==3
        indexSaveFileFolder=[saveRootFolder,'\Index\DataAlignment'];
        indexSaveFileName=[indexSaveFileFolder,'\Index',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),'.txt'];
        if exist(indexSaveFileName,'file')%���ھ�ֱ�Ӷ�
            Index_Alignment=importdata(indexSaveFileName);
        else
            Index_Alignment=ObtainDataAlignmentIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth);
        end
    end
    %----------------------------------------------------------------------
    %��ѡ��������
    %�ж��ļ��Ƿ����
    if method==2||method==3
        indexSaveFileFolder=[saveRootFolder,'\Index\DataElect'];
        maxIndexSaveFileName=[indexSaveFileFolder,'\Index_Max',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),'.txt'];
        minIndexSaveFileName=[indexSaveFileFolder,'\Index_Min',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),'.txt'];
        if exist(maxIndexSaveFileName,'file')%���ھ�ֱ�Ӷ�
            Index_Max=importdata(maxIndexSaveFileName);
            Index_Min=importdata(minIndexSaveFileName);
        else
            [Index_Max,Index_Min]=ObtainDataElectIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth);
        end
    end
    for i_snr=1:length(SNR(:,1));
        snr=10^(0.1*SNR(i_snr));
        %------------------------------------------------------------------
        % Calculate RMSE, some matrices
        T_rmse =zeros(num_Frame,num_Montecarlo);
%         T_mse = zeros(1,num_Frame);
        Loss = 0; % Track Loss flag
        C1 = gridWidth^2; C2 = 2*gridWidth^2; % Track Loss Threashood (lower,upper)
        success = 0; % flag of tracks which are successfully tracked
        %------------------------------------------------------------------
        for i_Montecarlo=1:num_Montecarlo
            %----------------------------------------------------------------------
            for i_Frame=1:length(Time(:,1));
                time=Time(i_Frame,1);
                Index_Target=GetTargetIndex(time,Coord_TargetBegin,Velocity,Coord_T,distanceUnit,num_T,num_R,num_Target,num_Channel);
                Coord_Target=Coord_TargetBegin+time*repmat(Velocity,num_Target,1);
                %%
                %�������ͳ����
                Measure=GetMeasure(snr,num_distanceUnit,Index_Target,num_Target,num_Channel);
                %%
                %----------------------------------------------------------
                %������
                if method==1||method==3
                    Result_Alignment=DetectionOfDataAlignment(Index_Alignment,Measure,num_XCoord,num_YCoord,gate_Alignment,Area,Temp_XCoord,Temp_YCoord,gridWidth,num_Channel,i_Frame,Coord_Target,flag_FigureDetect);
                    saveFileName_AlignmentResult=[saveFolder_AlignmentResult,'\',num2str(i_Frame,'%02d'),'.txt'];
                    save(saveFileName_AlignmentResult,'Result_Alignment','-ascii')
                end
                %----------------------------------------------------------
                %��ѡ�����
                if method==2||method==3
                    Result_Elect=Copy_of_DetectionOfDataElect(Index_Min,Index_Max,Measure,num_XCoord,num_YCoord,gate_Alignment,Area,Temp_XCoord,Temp_YCoord,gridWidth,num_Channel,num_Grid,Coord_Target,flag_FigureDetect);
                    saveFileName_ElectResult=[saveFolder_ElectResult,'\',num2str(i_Frame,'%02d'),'.txt'];
                    save(saveFileName_ElectResult,'Result_Elect','-ascii')
                end

            end
            %--------------------------------------------------------------
            %��ͼ
%             filename1='../Figure/Detection.gif';
%             
%             for i_Frame=1:length(Time(:,1));
%                 im=imread(['../Figure/Detection',num2str(i_Frame),'.bmp']);
%                 [imind,map2] = rgb2ind(im,256);
%                 if i_Frame==1
%                     imwrite(imind,map2,filename1,'gif', 'Loopcount',inf,'DelayTime',0.2)
%                 else
%                     imwrite(imind,map2,filename1,'gif','writeMode','append','delaytime',0.2);
%                 end
%             end
            %--------------------------------------------------------------
            %%
            %����
            if flag_Track==1
                Target_Location=repmat(Coord_TargetBegin,length(Time),1)+Time*Velocity;%�����õ�Ŀ��
                if method==1||method==3
                    [l,T_mse]=MyTrack(Area,Target_Location,Velocity,gridWidth,saveFolder_AlignmentResult,flag_FigureTrack);
                end
                if method==2||method==3
                    [l,T_mse]=MyTrack(Area,Target_Location,Velocity,gridWidth,saveFolder_ElectResult,flag_FigureTrack);
                end
                if l == 0
                    T_rmse(:,i_Montecarlo) = T_mse;%ÿ�з�һ�����ؿ����֡�����
                    % MIU(:,:,ii) = Miu; %���ף��ң�
                    success = success + 1;
                else
                    Loss = Loss + 1;
                end
            end
        end
        %% Statisticsͳ��
        %[clutter_x,clutter_y] = find(data_detection); % Radar clutter plots
        %���޸�1����˧-20150124-ԭ���������ƽ����ѭ���⣨�����ã��Ҹĵ���ͼ������
%         RMSE = zeros(1,num_Frame);
        % % Miu_aver = zeros(3,total_frame); %���ף�
        %Miu_aver = zeros(2,total_frame); %���ң�
        if flag_Track==1
            fprintf('GridWidth=%d,success=%d',GridWidth(i_GridWidth),success)
            for i_Frame = 1:num_Frame
                RMSE(i_GridWidth,i_Frame) = sqrt(sum(T_rmse(i_Frame,:))/success);
                % Miu_aver(1,ii) = sum(MIU(1,ii,:))/success; % ���ף��ң�
                %Miu_aver(2,ii) = sum(MIU(2,ii,:))/success; % ���ף��ң�
                % %     Miu_aver(3,ii) = sum(MIU(3,ii,:))/success; % ���ף�
            end
            figure(figure_RMSE)
            plot(RMSE(i_GridWidth,:),Mark_RMSE{i_GridWidth,:})
            hold on
            axis on
%             title('The Root Mean Square Error Curve');
            title('�������');
            xlabel('Frame/֡');
            ylabel('RMSE/m');
            Legend_RMSE{i_GridWidth,1}=['դ����=',num2str(gridWidth),'m'];
        end
    end
end
if flag_Track==1
    legend(Legend_RMSE)
end
fprintf('End time = %s\n',datestr(now))
