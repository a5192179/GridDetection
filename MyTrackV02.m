function [l,T_mse,Xpda,xfilt]=MyTrackV02(Area,Target_Location,Velocity,gridWidth,Result_Detection,flag_FigureTrack)
% global i_Figure
T=1;
total_frame = length(Target_Location(:,1));
T_mse = zeros(1,total_frame);

C1 = gridWidth^2; C2 = 2*gridWidth^2;
ss = 4; % state size
os = 2; % observation size
F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
H = [1 0 0 0; 0 1 0 0];
Q = 0.1*eye(ss);
R = gridWidth*eye(os);
initx = [Target_Location(1,:),Velocity]';
initV = 10*eye(ss);
One_Result_Detection=zeros(total_frame,2);
%把可能存在的多个目标简化成1个目标

for i_Frame = 1:total_frame
    Result_Detection_ThisFrame=zeros(1,2);
    Result_Detection_ThisFrame(1,1:2)=Result_Detection(i_Frame,1,:);
%     Result_Detection_ThisFrame=Result_Detection_ThisFrame;
    if isempty(Result_Detection_ThisFrame)
        Result_Detection_LastFrame=Result_Detection(i_Frame-1,:,:);
        Result_Detection_LastFrame(Result_Detection_LastFrame==1.211)=[];
        Result_Detection_ThisFrame=[Result_Detection_LastFrame(1,1)+Velocity(1,1)*T,Result_Detection_LastFrame(1,2)+Velocity(1,2)*T];
    end
    One_Result_Detection(i_Frame,:)=Result_Detection_ThisFrame(1,:);
end
[xfilt, Vfilt, VVfilt, loglik] = kalman_filter(One_Result_Detection', F, H, Q, R, initx, initV);
k = 0; l = 0;
for i_Frame = 1:total_frame
    T_mse(:,i_Frame) = (xfilt(1,i_Frame)-Target_Location(i_Frame,1))^2 + (xfilt(2,i_Frame)-Target_Location(i_Frame,2))^2; % square error per frame
    
    if sqrt(T_mse(:,i_Frame)) > C1 & sqrt(T_mse(:,i_Frame)) < C2%一般丢失
        k = k+1;
    elseif sqrt(T_mse(:,i_Frame))<C1
        k = 0;
    elseif sqrt(T_mse(:,i_Frame)) > C2%严重丢失，这次跟踪失败
        l = l + 1;
    end
    if k == 8
        k = 0; l = l + 1;
    end
end
Xpda=xfilt(:,end)';
%%
% if flag_FigureTrack==1
%     figure_TrackPoint=figure(i_Figure);
%     i_Figure=i_Figure+1;
%     figure(figure_TrackPoint)
%     for i_Frame=1:total_frame
%         plot_True=plot(Target_Location(:,1),Target_Location(:,2),'ro-');
%         %     legend(plot_True,'True')
%         hold on
%         %     if i>2
%         %         [clutter_x,clutter_y] = find(Save_data_detection(:,:,i));% Radar clutter plots
% %         detectionResultSaveFileName=[ResultSaveFolder,'\',num2str(i_Frame,'%02d'),'.txt'];
% %         Result_Detection=importdata(detectionResultSaveFileName);
% 
%         plot_Detect=plot(One_Result_Detection(:,1),One_Result_Detection(:,2),'m*');
%         %     legend(plot_Detect,'Detect')
%         hold on
%         %     end
%         
%         plot_Estimate=plot(xfilt(1,1:i_Frame),xfilt(2,1:i_Frame),'b.-');
%         legend([plot_True,plot_Detect,plot_Estimate],'True','Detect','Estimate','Location','NorthWest')
%         hold on
%         
%         %     hold off
%         title(['跟踪轨迹 栅格宽度=',num2str(gridWidth),'m']);
%         xlabel('X/m');
%         ylabel('Y/m');
%         grid on
%         axis([Area(1,1) Area(1,2) Area(2,1) Area(2,2)])
%         %     legend('Estimate','True');
% %         print(gcf,'-dbmp',strcat('../Figure/lsmycv',num2str(i_Frame),'.bmp'));
%     end
% end