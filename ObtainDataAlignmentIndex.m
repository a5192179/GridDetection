%**************************************************************************
%*��Creat time����2017-01-27 15:08          ��Version����0.0
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*�������ݶ�������
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
function Index=ObtainDataAlignmentIndex(Coord_T,distanceUnit,Temp_XCoord,Temp_YCoord,num_R,num_T,Coord_R,num_Channel,saveRootFolder,order_Radar,order_Aera,gridWidth)
% Temp_XCoord=(Area(1,1):gridWidth:Area(1,2))';
% num_XCoord=length(Temp_XCoord(:,1));
% Temp_YCoord=(Area(2,1):gridWidth:Area(2,2))';
% num_YCoord=length(Temp_YCoord(:,1));
%��չ��XY����һһ��Ӧ
Xcoord=kron(Temp_XCoord,ones(length(Temp_YCoord),1));
Ycoord=repmat(Temp_YCoord,length(Temp_XCoord),1);
Point_Gird=[Xcoord+0.5*gridWidth,Ycoord-0.5*gridWidth];
num_Grid=length(Point_Gird(:,1));
Index=zeros(num_Channel,num_Grid);
for i_T=1:num_T
    for i_R=1:num_R
        %����
        Temp_Radar_T=Coord_T(i_T,:);
        %��չ��һ�������״��Ӧ����դ��
        Point_Radar_T=repmat(Temp_Radar_T,num_Grid,1);
        %���㷢�����
        Distance_T=Distance2(Point_Radar_T,Point_Gird);
        %����
        Temp_Radar_R=Coord_R(i_R,:);
        %��չ��һ�������״��Ӧ����դ��
        Point_Radar_R=repmat(Temp_Radar_R,num_Grid,1);
        %������վ���
        Distance_R=Distance2(Point_Radar_R,Point_Gird);
        Distance=Distance_T+Distance_R;
        Order_Unit=(ceil(Distance/distanceUnit))';
        Index((i_T-1)*num_R+i_R,:)=(Order_Unit-1)*num_Channel+(i_T-1)*num_R+i_R;
%         Index((i_T-1)*num_R+i_R,:)=(ceil(Distance/distanceUnit))';
    end
end
if ~strcmp(saveRootFolder,'-1')
    saveFileFolder=[saveRootFolder,'\Index\DataAlignment'];
    if ~exist(saveFileFolder,'dir')
        mkdir(saveFileFolder)
    end
    saveFileName=[saveFileFolder,'\Index',num2str(order_Radar,'%02d'),num2str(gridWidth,'%03d'),num2str(order_Aera,'%02d'),'.txt'];%�����ţ�դ���С��ţ��״�������ţ��״�λ�ñ��
    save(saveFileName,'Index','-ascii')
end


