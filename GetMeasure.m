%**************************************************************************
%*【Creat time】：2017-01-28 14:35          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*产生回波信号
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
function Measure=GetMeasure(snr,num_distanceUnit,Index_Target,num_Target,num_Channel)
Measure=abs(sqrt(1/2)*(randn(num_Channel,num_distanceUnit)+1i*randn(num_Channel,num_distanceUnit))).^2;
for i_Channel=1:num_Channel
    for i_Target=1:num_Target
        Col=Index_Target((i_Channel-1)*num_Target+i_Target,1)-2:Index_Target((i_Channel-1)*num_Target+i_Target,1)+2;
        Measure(i_Channel,Col)=Measure(i_Channel,Col)+abs((sqrt(1/2)*(sqrt(10^(snr/10))+sqrt(10^(snr/10))*1i)).^2)*[0.3,0.6,1,0.6,0.3];
    end
end