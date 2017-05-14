%**************************************************************************
%*【Creat time】：2017-01-27 15:08          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*管理主函数的输入参数
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
function [Coord_T,Coord_R,GridWidth,Area]=SelectParameter(order_Radar,order_GridWidth,order_Aera)
%%
switch order_Radar
    case 1
        Coord_T=[
            0 50;%发射雷达坐标1 2 3 4 5 XY
            0 25;
            0 0;
            0 -25;
            0 -50]*1e3;
        Coord_R=Coord_T;
    case 2
        Coord_T=[
            0 50;%发射雷达坐标1 2 3 XY
            0 0;
            0 -25]*1e3;
        Coord_R=Coord_T;
    case 3
        Coord_T=[
            0 0;
            0 -25]*1e3;
        Coord_R=Coord_T;
    case 4
        Coord_T=[
            25 50;%发射雷达坐标1 2 3 4 5 XY
            12.5 25;
            0 0;
            -12.5 -25;
            -25 -50]*1e3;
        Coord_R=Coord_T;
    case 5
        Coord_T=[
            0 0;
            50 -50]*1e3;
        Coord_R=Coord_T;
    case 6
        Coord_T=[
            0 0;
            0 -25]*1e3;
        Coord_R=Coord_T(1,:);
    case 7
        Coord_T=[
            0 50;%发射雷达坐标1 2 3 4 5 XY
            0 25;
            0 0;
            0 -25;
            0 -50]*1e3;
        Coord_R=Coord_T(1,:);
    otherwise
        disp('Wrong radar order.')
        
end
%%
switch order_GridWidth
    case 1
        GridWidth=50;
    case 2
        GridWidth=10;
    case 3
        GridWidth=100;
    case 4
        GridWidth=[10;30;50];
    case 5
        GridWidth=30;
    case 6
        GridWidth=[50;100;200];%[10;20;50];%[50;100;200;300];
    case 7
        GridWidth=[10;30;50;70;90;110;130;150;170;200];
    case 8
%         GridWidth=[40;50;60;70;100;150;200];
%         GridWidth=[10;20;30;40;50;60;70;90;130;150;200];
        GridWidth=(10:5:150)';%[10;2;3;5;8;10;15;30;100;200];%[40;45;48;50;52;55;60];
%         GridWidth=[50;100;200];
    otherwise
        disp('Wrong gird width order.')
end
%%
switch order_Aera
    case 1
        Area=[49.5,50.5;          %X
            -0.5,0.5]*1e3;      %Y
    case 2
        Area=[47.5,52.5;          %X
            -2.5,2.5]*1e3;      %Y
    case 3
        Area=[49,51;          %X
            -1,1]*1e3;      %Y
    case 4
        Area=[45,55;          %X
            -5,5]*1e3;      %Y
    case 5
        Area=[22.5,27.5;          %X
            -12.5,-7.5]*1e3;      %Y
    case 6
        Area=[47.5,52.5;          %X
            -7.5,-2.5]*1e3;      %Y
    case 7
        Area=[49.5,50.5;          %X
            -1.5,-0.5]*1e3;      %Y
    case 8
        Area=[25.25,25.5;          %X
            -1.25,-1]*1e3;      %Y
    case 9
        Area=[47.5,52.5;          %X
            -3.5,1.5]*1e3;      %Y
    otherwise
        disp('Wrong area order.')
end
