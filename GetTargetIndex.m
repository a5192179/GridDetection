%**************************************************************************
%*【Creat time】：2017-01-28 14:35          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*产生目标在各通道的索引，用来给回波加信号，没有考虑目标在最大距离单元之外的情况
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
function Index_Target=GetTargetIndex(time,Coord_TargetBegin,Velocity,Coord_T,distanceUnit,num_T,num_R,num_Target,num_Channel)
Coord_R=Coord_T;
Point_Radar=[kron(Coord_T,ones(num_R,1));repmat(Coord_R,num_T,1)];%列排列的雷达坐标，前面一半为发射111222333，后面一半为接收123123123,通道11通道12通道13通道21通道22通道23
%计算目标位置
% Coord_Target=Coord_TargetBegin;
Coord_Target=Coord_TargetBegin+time*repmat(Velocity,num_Target,1);
%扩展雷达坐标，一个通道对应的所有目标
Point_Radar=kron(Point_Radar,ones(num_Target,1));%c1c1c1c2c2c2
Point_Target=repmat(repmat(Coord_Target,num_Channel,1),2,1);%t1t2t3t1t2t3,前一半为发射，后一半为接收
Distance=Distance2(Point_Radar,Point_Target);
Distance=Distance(1:end/2)+Distance(end/2+1:end);
Index_Target=ceil(Distance/distanceUnit);%通道1的所有目标索引；通道2的所有目标索引;...