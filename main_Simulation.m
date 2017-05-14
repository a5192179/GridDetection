%**************************************************************************
%*【Creat time】：2017-01-27 15:08          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*仿真程序的主函数，包含
%*
%*
%*
%*【Description】：
%*
%*
%*-------------------------------------------------------------------------
%*【Modification】：****-**-** **：**       【Version】：*.*
%*
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*
%*
%*
%*【Description】：
%*
%*
%*
%**************************************************************************
clc
close all
clear all
fprintf('Begin time = %s\n',datestr(now))
%%
%参数设置
saveRootFolder='..\Data';
saveFolder_AlignmentResult=[saveRootFolder,'\DetectionResult\Alignment'];
if ~exist(saveFolder_AlignmentResult,'dir') 
    mkdir(saveFolder_AlignmentResult)
end 
saveFolder_ElectResult=[saveRootFolder,'\DetectionResult\Elect'];
if ~exist(saveFolder_ElectResult,'dir') 
    mkdir(saveFolder_ElectResult)
end 
method=1;%1对齐，2优选，3都有
flag_Track=1;
flag_FigureDetect=0;
flag_FigureDetect_Gif=0;
flag_FigureTrack=0;%静态栅格画静态栅格的单帧实时图，跟踪栅格画跟踪栅格的单帧实时图
flag_Track_Gif=0;%画静态跟踪的动图
flag_TrackGrid=1;%启动跟踪栅格
flag_FigureTrackGridandArea=1;%统计跟踪栅格参数，在单帧实时图画上跟踪栅格图
flag_TrackGrid_Gif=0;%画跟踪栅格的动图
rhoA=8;
rhoG=1/10;
%雷达位置编号，栅格大小编号，区域编号
order_Radar=3;
order_GridWidth=6;
order_Aera=9;
[Coord_T,Coord_R,GridWidth,Area]=SelectParameter(order_Radar,order_GridWidth,order_Aera);
Area_Origin=Area;
%雷达参数
% Coord_T=[
%     0 50;%发射雷达坐标1 2 3 4 5 XY
%     0 25;
%     0 0;
%     0 -25;
%     0 -50]*1e3;
% Coord_R=Coord_T;
distanceUnit=200;%单位：m
detectRange=300e3;%雷达最远探测距离
num_distanceUnit=ceil(detectRange/distanceUnit);
%仿真场景
% Area=[49.5,50.5;          %X
%     -0.5,0.5]*1e3;      %Y
%
% Coord_TargetBegin=[50.4,0.4]*1e3;
Coord_TargetBegin=[52,2]*1e3;
% Coord_TargetBegin=[50.2,0.2]*1e3;
%
% Coord_target=Coord_targetBegin;
% Velocity=[-50,-50];
Time=(0:79)';%79
num_Frame=length(Time);
% Time=0;
% TimeBegin=zeros(size(Time));
% Coord_target(:,1)=Coord_target(:,1)+velocity*(Time-TimeBegin)
% GirdWidth=[30;40;50;60;70;100;150];%单位
% GridWidth=50;
%仿真参数
% SNR=[linspace(-10,-2,5)';linspace(0,15,31)';linspace(16,20,5)'];
SNR=14;
num_Montecarlo=1;
%%
%计算数量
num_T=length(Coord_T(:,1));
num_Target=length(Coord_TargetBegin(:,1));
num_R=length(Coord_R(:,1));
num_Channel=num_T*num_R;
num_GridWidth=length(GridWidth(:,1));
num_TargetMax=100;
F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
% global i_Figure
i_Figure=1;
%%
%载入门限
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
%产生索引
%%
%--------------------------------------------------------------------------
%跟踪画图参数
if flag_Track==1
    figure_RMSE=figure(i_Figure);
    i_Figure=i_Figure+1;
    Mark_RMSE={'b-';'r:';'m-.';'k--';'g-'};
    Legend_RMSE=cell(num_GridWidth,1);
    if flag_FigureTrackGridandArea==1
        Legend_GridWidth=cell(num_GridWidth,1);
    end
    RMSE = zeros(length(GridWidth(:,1)),num_Frame);
end
%--------------------------------------------------------------------------

for i_GridWidth=1:num_GridWidth

    for i_SNR=1:length(SNR(:,1));
        snr=SNR(i_SNR);
        %------------------------------------------------------------------
        % Calculate RMSE, some matrices
        T_rmse =zeros(num_Frame,num_Montecarlo);
%         T_mse = zeros(1,num_Frame);
        Loss = 0; % Track Loss flag
        
        success = 0; % flag of tracks which are successfully tracked
        %------------------------------------------------------------------
        for i_Montecarlo=1:num_Montecarlo
            %重置Area
            Area=Area_Origin;
            %重置saveRootFolder
            saveRootFolder='..\Data';
            %设置目标
            if (flag_FigureDetect||flag_FigureTrack||flag_FigureTrackGridandArea||flag_FigureDetect)&&i_Montecarlo==1
                theta=45;
            else
                theta=360*rand(1);
            end
            r=(Area(1,2)-Area(1,1))/2;
            x_TargetBegin=r*cosd(theta)+(Area(1,2)+Area(1,1))/2;
            y_TargetBegin=r*sind(theta)+(Area(2,2)+Area(2,1))/2;
            x_TargetEnd=Area(1,2)+Area(1,1)-x_TargetBegin;
            y_TargetEnd=Area(2,2)+Area(2,1)-y_TargetBegin;
            Coord_TargetBegin=[x_TargetBegin,y_TargetBegin];
            Velocity=[(x_TargetEnd-x_TargetBegin)/num_Frame,(y_TargetEnd-y_TargetBegin)/num_Frame];
            %--------------------------------------------------------------
            gridWidth=GridWidth(i_GridWidth,1);%单位
            C1 = gridWidth^2; C2 = 2*gridWidth^2; % Track Loss Threashood (lower,upper)
            %划分栅格
            Temp_XCoord=(Area(1,1):gridWidth:Area(1,2))';%从左到右
            num_XCoord=length(Temp_XCoord(:,1));
            Temp_YCoord=(Area(2,2):-gridWidth:Area(2,1))';%从上到下
            num_YCoord=length(Temp_YCoord(:,1));
            num_Grid=num_XCoord*num_YCoord;
            %-----------------------------------------------------------
            %对齐载入索引
            %判断文件是否存在
            if method==1||method==3
                indexSaveFileFolder=[saveRootFolder,'\Index\DataAlignment'];
                indexSaveFileName=[indexSaveFileFolder,'\Index',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.txt'];
                if exist(indexSaveFileName,'file')&&~(strcmp(num2str(gridWidth,'%03d'),'010')&&strcmp(num2str(order_Aera,'%02d'),'04'))%存在就直接读
                    Index_Alignment=importdata(indexSaveFileName);
                else
                    Index_Alignment=ObtainDataAlignmentIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth);
                end
            end
            %-----------------------------------------------------------
            %优选载入索引
            %判断文件是否存在
            if method==2||method==3
                indexSaveFileFolder=[saveRootFolder,'\Index\DataElect'];
                maxIndexSaveFileName=[indexSaveFileFolder,'\Index_Max',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.txt'];
                minIndexSaveFileName=[indexSaveFileFolder,'\Index_Min',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.txt'];
                if exist(maxIndexSaveFileName,'file')%存在就直接读
                    Index_Max=importdata(maxIndexSaveFileName);
                    Index_Min=importdata(minIndexSaveFileName);
                else
                    [Index_Max,Index_Min]=ObtainDataElectIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth);
                end
            end
            %--------------------------------------------------------------
            Result_Alignment=zeros(num_Frame,num_TargetMax,2)+1.211;
            Result_Elect=zeros(num_Frame,num_TargetMax,2)+1.211;
            %----------------------------------------------------------------------
            %初始化跟踪栅格的画图的参数
            if i_Montecarlo==1&&flag_FigureTrackGridandArea==1
                Area_AllFrame=zeros(num_Frame,5);%每行四个顶点，1个栅格宽度
                Area_AllFrame(1,:)=[Area(1,1) Area(1,2) Area(2,1) Area(2,2) gridWidth];%每行四个顶点，1个栅格宽度
                Area_AllFrame(2,:)=[Area(1,1) Area(1,2) Area(2,1) Area(2,2) gridWidth];%每行四个顶点，1个栅格宽度
            end
            %--------------------------------------------------------------
            for i_Frame=1:length(Time(:,1));
                time=Time(i_Frame,1);
                Index_Target=GetTargetIndex(time,Coord_TargetBegin,Velocity,Coord_T,distanceUnit,num_T,num_R,num_Target,num_Channel);
                Coord_Target=Coord_TargetBegin+time*repmat(Velocity,num_Target,1);
                %%
                %产生检测统计量
                Measure=GetMeasure(snr,num_distanceUnit,Index_Target,num_Target,num_Channel);
                %%
                %----------------------------------------------------------
                %对齐检测
                if method==1
                    temp_Result_Alignment=DetectionOfDataAlignment(Index_Alignment,Measure,num_XCoord,num_YCoord,gate_Alignment,Area,Temp_XCoord,Temp_YCoord,gridWidth,num_Channel,i_Frame,Coord_Target,flag_FigureDetect,flag_FigureDetect_Gif,i_Montecarlo);
                    num_Result=length(temp_Result_Alignment(:,1));
                    if num_Result<=num_TargetMax
                        Result_Alignment(i_Frame,1:num_Result,:)=temp_Result_Alignment;
                    else
                        Result_Alignment(i_Frame,:,:)=temp_Result_Alignment(1:100,:);%最大100个目标
                    end
%                     saveFileName_AlignmentResult=[saveFolder_AlignmentResult,'\',num2str(i_Frame,'%02d'),'.txt'];
%                     save(saveFileName_AlignmentResult,'Result_Alignment','-ascii')
                end
                %----------------------------------------------------------
                %优选检测检测
                if method==2
                    temp_Result_Elect=DetectionOfDataElect(Index_Min,Index_Max,Measure,num_XCoord,num_YCoord,gate_Alignment,Area,Temp_XCoord,Temp_YCoord,gridWidth,num_Channel,num_Grid,Coord_Target,flag_FigureDetect,flag_FigureDetect_Gif,i_Montecarlo,i_Frame);
                    num_Result=length(temp_Result_Elect(:,1));
                    if num_Result<=num_TargetMax
                        Result_Elect(i_Frame,1:num_Result,:)=temp_Result_Elect;
                    else
                        Result_Elect(i_Frame,:,:)=temp_Result_Elect(1:100,:);%最大100个目标
                    end
%                     saveFileName_ElectResult=[saveFolder_ElectResult,'\',num2str(i_Frame,'%02d'),'.txt'];
%                     save(saveFileName_ElectResult,'Result_Elect','-ascii')
                end
                %----------------------------------------------------------
                %画检测结果动态图
                if i_Montecarlo==1&&flag_FigureDetect_Gif==1&&i_Frame==num_Frame
                    if method==1
                        filename1=['../Figure/Detect/Alignment',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.gif'];
                    else
                        filename1=['../Figure/Detect/Elect',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.gif'];
                    end
                    
                    for temp_I_Frame=1:length(Time(:,1));
                        im=imread(['../Figure/Detect/Detect',num2str(temp_I_Frame),'.bmp']);
                        [imind,map2] = rgb2ind(im,256);
                        if temp_I_Frame==1
                            imwrite(imind,map2,filename1,'gif', 'Loopcount',inf,'DelayTime',0.2)
                        else
                            imwrite(imind,map2,filename1,'gif','writeMode','append','delaytime',0.2);
                        end
                    end
                    print(gcf,'-dbmp',strcat('../Figure/lsmycv',num2str(temp_I_Frame),'.bmp'));
                end
                %----------------------------------------------------------
                %跟踪栅格跟踪
                if flag_TrackGrid==1&&i_Frame>=3
                    Target_Location=repmat(Coord_TargetBegin,i_Frame,1)+Time(1:i_Frame)*Velocity;%仅适用单目标
                    if method==1
%                         [l,T_mse,Xk]=MyTrack(Area,Target_Location,Velocity,gridWidth,Result_Alignment(1:i_Frame,:,:),flag_FigureTrack);
                        [l,T_mse,Xpda,xfilt]=MyTrackV02(Area,Target_Location,Velocity,gridWidth,Result_Alignment(1:i_Frame,:,:),flag_FigureTrack);
                    end
                    if method==2
%                         [l,T_mse,Xk]=MyTrack(Area,Target_Location,Velocity,gridWidth,Result_Elect(1:i_Frame,:,:),flag_FigureTrack);
                        [l,T_mse,Xpda,xfilt]=MyTrackV02(Area,Target_Location,Velocity,gridWidth,Result_Elect(1:i_Frame,:,:),flag_FigureTrack);
                    end
                    %--------------------------------------------------
                    %每个栅格宽度第一次模拟最后一帧开新图画
                    if method==1
                        Result_Detect=Result_Alignment;
                    else
                        Result_Detect=Result_Elect;
                    end
                    if i_Montecarlo==1&&i_Frame==num_Frame&&flag_FigureTrack==1
                        for temp_I_Frame=1:num_Frame
                            if temp_I_Frame==1
                                figure_TrackPoint=figure(i_Figure);
                                i_Figure=i_Figure+1;
                                figure(figure_TrackPoint)
                            end
                            %画真实点迹
                            plot_True=plot(Target_Location(1:temp_I_Frame,1),Target_Location(1:temp_I_Frame,2),'ro-');
                            hold on
                            %------------------------------------------
                            %找出每一帧的目标，画上去
                            Result_Detection_ThisFrame=zeros(100,2);
                            Result_Detection_ThisFrame(:,1:2)=Result_Detect(temp_I_Frame,:,:);
                            Result_Detection_ThisFrame(Result_Detection_ThisFrame==1.211)=[];
                            if ~isempty(Result_Detection_ThisFrame)
                                plot_Detect=plot(Result_Detection_ThisFrame(:,1),Result_Detection_ThisFrame(:,2),'m*');
                                hold on
                            end
                            %------------------------------------------
                            %预测结果
                            plot_Estimate=plot(xfilt(1,1:temp_I_Frame),xfilt(2,1:temp_I_Frame),'b.-');
                            hold on
                            if flag_FigureTrackGridandArea==1
                                LinePoint_Area=[Area_AllFrame(temp_I_Frame,1) Area_AllFrame(temp_I_Frame,3);
                                    Area_AllFrame(temp_I_Frame,2) Area_AllFrame(temp_I_Frame,3);
                                    Area_AllFrame(temp_I_Frame,2) Area_AllFrame(temp_I_Frame,4);
                                    Area_AllFrame(temp_I_Frame,1) Area_AllFrame(temp_I_Frame,4);
                                    Area_AllFrame(temp_I_Frame,1) Area_AllFrame(temp_I_Frame,3)];
                                plot_Area=plot(LinePoint_Area(:,1),LinePoint_Area(:,2),'g-');
                                hold on
                                legend([plot_True,plot_Detect,plot_Estimate,plot_Area],'真实目标','检测结果','跟踪轨迹','检测区域','Location','NorthWest')
                            else%不画跟踪栅格的区域
                                legend([plot_True,plot_Detect,plot_Estimate],'真实目标','检测结果','跟踪轨迹','Location','NorthWest')
                            end
                            
                            %     hold off
                            title(['动态栅格跟踪轨迹 栅格宽度=',num2str(GridWidth(i_GridWidth)),'m']);
                            xlabel('X/m');
                            ylabel('Y/m');
                            grid on
                            axis([Area_Origin(1,1) Area_Origin(1,2) Area_Origin(2,1) Area_Origin(2,2)])
                            if flag_TrackGrid_Gif==1
                                print(gcf,'-dbmp',strcat('../Figure/Track/Track',num2str(temp_I_Frame),'.bmp'));
                            end
                        end
                        %画动图
                        if flag_TrackGrid_Gif==1
                            filename1=['../Figure/Track/DynamicGrid',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.gif'];
                            
                            for temp_I_Frame=1:length(Time(:,1));
                                im=imread(['../Figure/Track/Track',num2str(temp_I_Frame),'.bmp']);
                                [imind,map2] = rgb2ind(im,256);
                                if temp_I_Frame==1
                                    imwrite(imind,map2,filename1,'gif', 'Loopcount',inf,'DelayTime',0.2)
                                else
                                    imwrite(imind,map2,filename1,'gif','writeMode','append','delaytime',0.2);
                                end
                            end
                            print(gcf,'-dbmp',strcat('../Figure/lsmycv',num2str(temp_I_Frame),'.bmp'));
                        end
                    end
                    %画栅格尺寸变化
                    if i_Montecarlo==1&&i_Frame==num_Frame&&flag_FigureTrackGridandArea==1
                        if i_GridWidth==1
                            figure_GridWidth=figure(i_Figure);
                            i_Figure=i_Figure+1;
                        end
                        figure(figure_GridWidth)
                        plot_True=plot(Time(1:16,1),Area_AllFrame(1:16,5),Mark_RMSE{i_GridWidth,:},'LineWidth',2);
                        axis([0 15 0 210])
                        hold on
                        title('栅格宽度变化');
                        xlabel('帧');
                        ylabel('栅格宽度/m');
                        grid on
                        Legend_GridWidth{i_GridWidth,1}=['栅格宽度=',num2str(GridWidth(i_GridWidth)),'m'];
                        if i_GridWidth==num_GridWidth
                            legend(Legend_GridWidth)
                        end
                    end
                    %画区域尺寸变化
                    if i_Montecarlo==1&&i_Frame==num_Frame&&flag_FigureTrackGridandArea==1
                        if i_GridWidth==1
                            figure_AreaWidth=figure(i_Figure);
                            i_Figure=i_Figure+1;
                        end
                        figure(figure_AreaWidth)
                        plot_True=plot(Time(1:16,1),Area_AllFrame(1:16,2)-Area_AllFrame(1:16,1),Mark_RMSE{i_GridWidth,:},'LineWidth',2);
                        hold on
                        title('检测区域宽度变化');
                        xlabel('帧');
                        ylabel('区域宽度/m');
                        grid on
                        Legend_GridWidth{i_GridWidth,1}=['栅格宽度=',num2str(GridWidth(i_GridWidth)),'m'];
                        if i_GridWidth==num_GridWidth
                            legend(Legend_GridWidth)
                        end
                    end
                    %动态栅格画图结束
                    %--------------------------------------------------
                    Xkadd1=F*Xpda'; % 一步预测状态向量,x位置，x速度，y位置，y速度
                    %------------------------------------------------------
                    %更新Area
                    Center_Area=Xkadd1([1,2]);
                    length_Area=1.2*sqrt(Xkadd1(3)^2+Xkadd1(4)^2);
                    if length_Area<gridWidth*rhoA
                        length_Area=gridWidth*rhoA;
                    end
                    Area=[Center_Area(1)-length_Area,Center_Area(1)+length_Area;Center_Area(2)-length_Area,Center_Area(2)+length_Area];
                    %更新栅格宽度
                    gridWidth=length_Area*rhoG;%单位
                    if gridWidth<50%distanceUnit/num_Channel
                        gridWidth=50;%distanceUnit/num_Channel;
                    elseif gridWidth>distanceUnit
                        gridWidth=distanceUnit;  
                    end
                    %------------------------------------------------------
                    %跟踪栅格的画图的参数
                    if i_Montecarlo==1&&flag_FigureTrackGridandArea==1
                        Area_AllFrame(i_Frame,:)=[Area(1,1) Area(1,2) Area(2,1) Area(2,2) gridWidth];%每行四个顶点，1个栅格宽度
                    end
                    %------------------------------------------------------
                    %划分栅格
                    Temp_XCoord=(Area(1,1):gridWidth:Area(1,2))';%从左到右
                    num_XCoord=length(Temp_XCoord(:,1));
                    Temp_YCoord=(Area(2,2):-gridWidth:Area(2,1))';%从上到下
                    num_YCoord=length(Temp_YCoord(:,1));
                    num_Grid=num_XCoord*num_YCoord;
                    %---------------------------------------------------
                    %对齐载入索引
                    saveRootFolder='-1';
                    %判断文件是否存在
                    
                    if method==1||method==3
                        Index_Alignment=ObtainDataAlignmentIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth);
                    end
                    %---------------------------------------------------
                    %优选载入索引
                    %判断文件是否存在
                    if method==2||method==3
                        [Index_Max,Index_Min]=ObtainDataElectIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth);
                    end
                end
                %----------------------------------------------------------

            end
            %--------------------------------------------------------------
            %动图
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
            %跟踪
            if flag_Track==1
                if flag_TrackGrid==0%前面没有跟踪
                    Target_Location=repmat(Coord_TargetBegin,i_Frame,1)+Time(1:i_Frame)*Velocity;%仅适用单目标
                    if method==1
                        [l,T_mse,Xpda,xfilt]=MyTrackV02(Area,Target_Location,Velocity,gridWidth,Result_Alignment,flag_FigureTrack);
%                         [l,T_mse,Xk]=MyTrack(Area,Target_Location,Velocity,gridWidth,Result_Alignment,flag_FigureTrack);
                    end
                    if method==2
                        [l,T_mse,Xpda,xfilt]=MyTrackV02(Area,Target_Location,Velocity,gridWidth,Result_Elect,flag_FigureTrack);
%                         [l,T_mse,Xk]=MyTrack(Area,Target_Location,Velocity,gridWidth,Result_Elect,flag_FigureTrack);
                    end
                    %------------------------------------------------------
                    %每个栅格宽度第一次模拟最后一帧开新图画，静态栅格实时跟踪画图
                    if method==1
                        Result_Detect=Result_Alignment;
                    else
                        Result_Detect=Result_Elect;
                    end
                    if i_Montecarlo==1&&i_Frame==num_Frame&&flag_FigureTrack==1
                        for temp_I_Frame=1:num_Frame
                            if temp_I_Frame==1
                                figure_TrackPoint=figure(i_Figure);
                                i_Figure=i_Figure+1;
                                figure(figure_TrackPoint)
                            end
                            plot_True=plot(Target_Location(1:temp_I_Frame,1),Target_Location(1:temp_I_Frame,2),'ro-');
                            hold on
                            %------------------------------------------
                            %找出每一帧的目标，画上去
                            Result_Detection_ThisFrame=zeros(100,2);
                            Result_Detection_ThisFrame(:,1:2)=Result_Detect(temp_I_Frame,:,:);
                            Result_Detection_ThisFrame(Result_Detection_ThisFrame==1.211)=[];
                            if ~isempty(Result_Detection_ThisFrame)
                                plot_Detect=plot(Result_Detection_ThisFrame(:,1),Result_Detection_ThisFrame(:,2),'m*');
                                hold on
                            end
                            %------------------------------------------
                            plot_Estimate=plot(xfilt(1,1:temp_I_Frame),xfilt(2,1:temp_I_Frame),'b.-');
                            legend([plot_True,plot_Detect,plot_Estimate],'真实目标','检测结果','跟踪轨迹','Location','NorthWest')
                            hold on
                            %     hold off
                            title(['静态栅格跟踪轨迹 栅格宽度=',num2str(gridWidth),'m']);
                            xlabel('X/m');
                            ylabel('Y/m');
                            grid on
                            axis([Area_Origin(1,1) Area_Origin(1,2) Area_Origin(2,1) Area_Origin(2,2)])
                            %     legend('Estimate','True');
                            if flag_Track_Gif==1
                                print(gcf,'-dbmp',strcat('../Figure/Track/Track',num2str(temp_I_Frame),'.bmp'));
                            end
                        end
                        %--------------------------------------------------
                        %制做动图
                        if flag_Track_Gif==1
                            filename1=['../Figure/Track/StaticGrid',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),num2str(distanceUnit,'%03d'),'.gif'];
                            for temp_I_Frame=1:length(Time(:,1));
                                im=imread(['../Figure/Track/Track',num2str(temp_I_Frame),'.bmp']);
                                [imind,map2] = rgb2ind(im,256);
                                if temp_I_Frame==1
                                    imwrite(imind,map2,filename1,'gif', 'Loopcount',inf,'DelayTime',0.2)
                                else
                                    imwrite(imind,map2,filename1,'gif','writeMode','append','delaytime',0.2);
                                end
                            end
                            print(gcf,'-dbmp',strcat('../Figure/lsmycv',num2str(temp_I_Frame),'.bmp'));
                        end
                        %--------------------------------------------------
                    end
                    %静态栅格实时跟踪画图结束
                    %------------------------------------------------------
                end
                if l == 0
                    T_rmse(:,i_Montecarlo) = T_mse;%每列放一次蒙特卡洛各帧的误差
                    % MIU(:,:,ii) = Miu; %（甲，乙）
                    success = success + 1;
                else
                    Loss = Loss + 1;
                end
            end
            if mod(i_Montecarlo,1)==0
                fprintf('SNR=%d,grid width=%d,i_Montecarlo=%d,time=%s\n',snr,gridWidth,i_Montecarlo,datestr(now))
            end
        end
        %% Statistics统计
        %[clutter_x,clutter_y] = find(data_detection); % Radar clutter plots
        %【修改1】李帅-20150124-原程序的噪声平面在循环外（这里）求得，我改到画图里面了
%         RMSE = zeros(1,num_Frame);
        % % Miu_aver = zeros(3,total_frame); %（甲）
        %Miu_aver = zeros(2,total_frame); %（乙）
        if flag_Track==1
            fprintf('GridWidth=%d,success=%d\n',GridWidth(i_GridWidth),success)
            for i_Frame = 1:num_Frame
                RMSE(i_GridWidth,i_Frame) = sqrt(sum(T_rmse(i_Frame,:))/success);
                % Miu_aver(1,ii) = sum(MIU(1,ii,:))/success; % （甲，乙）
                %Miu_aver(2,ii) = sum(MIU(2,ii,:))/success; % （甲，乙）
                % %     Miu_aver(3,ii) = sum(MIU(3,ii,:))/success; % （甲）
            end
            figure(figure_RMSE)
            plot(RMSE(i_GridWidth,:),Mark_RMSE{i_GridWidth,:},'LineWidth',2)
            hold on
            grid on
            axis on
%             title('The Root Mean Square Error Curve');
            title('跟踪误差');
            xlabel('帧');
            ylabel('RMSE/m');
            Legend_RMSE{i_GridWidth,1}=['栅格宽度=',num2str(GridWidth(i_GridWidth)),'m'];
        end
    end
end
if flag_Track==1
    legend(Legend_RMSE)
end
fprintf('End time = %s\n',datestr(now))
fprintf('End time = %s\n',datestr(now))