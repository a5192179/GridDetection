%**************************************************************************
%*��Creat time����2017-01-28 14:35          ��Version����0.0
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*����Ŀ���ڸ�ͨ�����������������ز����źţ�û�п���Ŀ���������뵥Ԫ֮������
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
function Index_Target=GetTargetIndex(time,Coord_TargetBegin,Velocity,Coord_T,distanceUnit,num_T,num_R,num_Target,num_Channel)
Coord_R=Coord_T;
Point_Radar=[kron(Coord_T,ones(num_R,1));repmat(Coord_R,num_T,1)];%�����е��״����꣬ǰ��һ��Ϊ����111222333������һ��Ϊ����123123123,ͨ��11ͨ��12ͨ��13ͨ��21ͨ��22ͨ��23
%����Ŀ��λ��
% Coord_Target=Coord_TargetBegin;
Coord_Target=Coord_TargetBegin+time*repmat(Velocity,num_Target,1);
%��չ�״����꣬һ��ͨ����Ӧ������Ŀ��
Point_Radar=kron(Point_Radar,ones(num_Target,1));%c1c1c1c2c2c2
Point_Target=repmat(repmat(Coord_Target,num_Channel,1),2,1);%t1t2t3t1t2t3,ǰһ��Ϊ���䣬��һ��Ϊ����
Distance=Distance2(Point_Radar,Point_Target);
Distance=Distance(1:end/2)+Distance(end/2+1:end);
Index_Target=ceil(Distance/distanceUnit);%ͨ��1������Ŀ��������ͨ��2������Ŀ������;...