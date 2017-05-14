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
function [Coord_T,GridWidth,Area]=SelectParameter(order_Radar,order_GridWidth,order_Aera)
%%
switch order_Radar
    case 1
        Coord_T=[
            0 50;%发射雷达坐标1 2 3 4 5 XY
            0 25;
            0 0;
            0 -25;
            0 -50]*1e3;
    case 2
        Coord_T=[
            0 50;%发射雷达坐标1 2 3 XY
            0 0;
            0 -25]*1e3;        
    case 3
        Coord_T=[
            0 0;
            0 -25]*1e3;
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
    otherwise
        disp('Wrong gird width order.')
end
%%
switch order_Aera
    case 1
        Area=[49.5,50.5;          %X
            -0.5,0.5]*1e3;      %Y
    otherwise
        disp('Wrong area order.')
end
