%**************************************************************************
%*��Creat time����2017-01-25 10:00          ��Version����0.0
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*
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
Coord_t=[
    0 50;%�����״�����1 2 3 4 5 XY
    0 25;
    0 0;
    0 -25;
    0 -50]*1e3;
Area=[49.5,50.5;          %X
      -0.5,0.5]*1e3;      %Y
BeginCoord=kron(Coord_t,ones(4,1));
EndCoord=repmat([Area(1,1),Area(2,1);Area(1,1),Area(2,2);Area(1,2),Area(2,1);Area(1,2),Area(2,2)],5,1);
EndCoord=repmat(mean(Area,2)',5,1);
Angle=AngleToEast(Coord_t,EndCoord);