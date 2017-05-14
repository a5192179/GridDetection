%**************************************************************************
%*【Creat time】：2017-02-06 17:02          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*生成数据对齐索引
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
function [Index_Max,Index_Min]=ObtainDataElectIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth)
% Temp_XCoord=(Area(1,1):gridWidth:Area(1,2))';
% num_XCoord=length(Temp_XCoord(:,1));
% Temp_YCoord=(Area(2,1):gridWidth:Area(2,2))';
% num_YCoord=length(Temp_YCoord(:,1));
%扩展到XY坐标一一对应
Xcoord=kron(Temp_XCoord,ones(length(Temp_YCoord),1));
Ycoord=repmat(Temp_YCoord,length(Temp_XCoord),1);
Point_Grid=[Xcoord,Ycoord;Xcoord,Ycoord-gridWidth;Xcoord+gridWidth,Ycoord;Xcoord+gridWidth,Ycoord-gridWidth];%左上，左下，右上，右下
num_Grid=length(Point_Grid(:,1))/4;
Index_Min=zeros(num_Channel,num_Grid);
Index_Max=zeros(num_Channel,num_Grid);
for i_T=1:num_T
    for i_R=1:num_R
        %发射
        Temp_Radar_T=Coord_T(i_T,:);
        %扩展到一个发射雷达对应所有栅格
        Point_Radar_T=repmat(Temp_Radar_T,4*num_Grid,1);
        %计算发射距离
        Distance_T=Distance2(Point_Radar_T,Point_Grid);
        %接收
        Temp_Radar_R=Coord_R(i_R,:);
        %扩展到一个接收雷达对应所有栅格
        Point_Radar_R=repmat(Temp_Radar_R,4*num_Grid,1);
        %计算接收距离
        Distance_R=Distance2(Point_Radar_R,Point_Grid);
        Distance=Distance_T+Distance_R;
        Order_Unit=(ceil(Distance/distanceUnit));
        Order_Unit=reshape(Order_Unit,num_Grid,4);
        Order_Unit_Max=max(Order_Unit,[],2);
        Order_Unit_Min=min(Order_Unit,[],2);
        Index_Max((i_T-1)*num_R+i_R,:)=(Order_Unit_Max-1)*num_Channel+(i_T-1)*num_R+i_R;
        Index_Min((i_T-1)*num_R+i_R,:)=(Order_Unit_Min-1)*num_Channel+(i_T-1)*num_R+i_R;
%         Index((i_T-1)*num_R+i_R,:)=(ceil(Distance/distanceUnit))';
    end
end
if ~strcmp(saveRootFolder,'-1')
    saveFileFolder=[saveRootFolder,'\Index\DataElect'];
    if ~exist(saveFileFolder,'dir')
        mkdir(saveFileFolder)
    end
    saveFileName=[saveFileFolder,'\Index_Max',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),'.txt'];%区域编号，栅格大小编号，雷达数量编号，雷达位置编号
    save(saveFileName,'Index_Max','-ascii')
    saveFileName=[saveFileFolder,'\Index_Min',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),'.txt'];%区域编号，栅格大小编号，雷达数量编号，雷达位置编号
    save(saveFileName,'Index_Min','-ascii')
end


