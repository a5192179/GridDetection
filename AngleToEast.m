%**************************************************************************
%*��Creat time����2016-05-21 19:18          ��Version����0.0
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*�������������������нǣ���ʱ�룬��ά�ѿ������꣬[0��,360�㣩
%*���룺
%*������Ƕ� angle ��λ��
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
% EndCoord=[1,1;-1,1;sqrt(3)/2,0.5;-sqrt(3)/2,0.5;-sqrt(3)/2,-0.5;sqrt(3)/2,-0.5];
% BeginCoord=zeros(size(EndCoord));
function Angle=AngleToEast(BeginCoord,EndCoord)
Distance=distance2(BeginCoord,EndCoord);
Asin=asind((EndCoord(:,2)-BeginCoord(:,2))./Distance);
CosAngle=(EndCoord(:,1)-BeginCoord(:,1))./Distance;%ԭ����acos_r=asind((Coord_r(i_r,2)-Coord_grid(1,2))/distance_r);
Angle=Asin.*(Asin>=0&CosAngle>0)+(180-Asin).*(Asin>0&CosAngle<=0)+(180-Asin).*(Asin<=0&CosAngle<0)+(2*180+Asin).*(Asin<0&CosAngle>=0);

function Distance=distance2(Point1,Point2)
Distance=sqrt((Point1(:,1)-Point2(:,1)).^2+(Point1(:,2)-Point2(:,2)).^2);