%**************************************************************************
%*¡¾Creat time¡¿£º2017-01-25 10:00          ¡¾Version¡¿£º0.0
%*¡¾Writer¡¿£ºLiShuai 461837929@qq.com
%*¡¾Writer department¡¿£º
%*¡¾Function¡¿£º
%*
%*
%*
%*
%*¡¾Description¡¿£º
%*
%*
%*-------------------------------------------------------------------------
%*¡¾Modification¡¿£º****-**-** **£º**       ¡¾Version¡¿£º*.*
%*
%*¡¾Writer¡¿£ºLiShuai 461837929@qq.com
%*¡¾Writer department¡¿£º
%*¡¾Function¡¿£º
%*
%*
%*
%*¡¾Description¡¿£º
%*
%*
%*
%**************************************************************************
Coord_t=[
    0 50;%·¢ÉäÀ×´ï×ø±ê1 2 3 4 5 XY
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