%**************************************************************************
%*【Creat time】：2017-01-29 21:56          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*检测，去鬼影
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
function Result=DetectionOfDataElect(Index_Min,Index_Max,Measure,num_XCoord,num_YCoord,gate_Alignment,Area,Temp_XCoord,Temp_YCoord,gridWidth,num_Channel,num_Grid,Coord_Target,flag_FigureDetect,flag_FigureDetect_Gif,i_Montecarlo,i_Frame)
gate_Elect=1.5*gate_Alignment;
num_Erasure=2;%前后各擦除2个单元
Detection=zeros(num_Channel,num_Grid);
%提取检测统计量
for i_Grid=1:num_Grid
    for i_Channel=1:num_Channel
        Detection(i_Channel,i_Grid)=max(Measure(Index_Min(i_Channel,i_Grid):num_Channel:Index_Max(i_Channel,i_Grid)));
    end
end
Detection=sum(Detection);
%检测
OriginalOrder=(1:num_XCoord*num_YCoord);%行向量，代表每个栅格的原始编号
Order_Result=Detection>gate_Elect;%每个过门限栅格在Detection中的编号
num_Result=sum(Order_Result);
Result=zeros(num_Result,2);
%--------------------------------------------------------------------------
%去鬼影前的点
Detection_Figure=reshape(Detection,num_YCoord,num_XCoord);
[Row,Col]=find(Detection_Figure>gate_Elect);
Coord_Result_X=Area(1,1)-0.5*gridWidth+Col*gridWidth;
Coord_Result_Y=Area(2,2)+0.5*gridWidth-Row*gridWidth;
% Temp_YCoord_Figure=zeros(size(Temp_YCoord));
% Detection_Figure=zeros(size(Detection));
% for i_Row=1:length(Temp_YCoord_Figure(:,1))
%     Temp_YCoord_Figure(i_Row,:)=Temp_YCoord(end-i_Row+1,:);
%     Detection_Figure(i_Row,:)=Detection(end-i_Row+1,:);
% end
%------------------------------------------------------------------------
if i_Montecarlo==1&&flag_FigureDetect==1
    temp_Figure_Detect=figure();
    imagesc(Temp_XCoord+gridWidth/2,Temp_YCoord-gridWidth/2,Detection_Figure)
    hold on
end
% % plot(50000,200,'y*','MarkerSize',10,'LineWidth',2)
% plot(Coord_Result_X,Coord_Result_Y,'g*','MarkerSize',10,'LineWidth',0.5)
%------------------------------------------------------------------------
%--------------------------------------------------------------------------
%去鬼影

for i_Result=1:num_Result
    %检测
    Order_Result=Detection>gate_Elect;%检测统计量里大于门限的位置
    if sum(Order_Result)==0
        i_Result=i_Result-1;
        break
    end

    %找到最大的
    [detection_Max,order_Max]=max(Detection);
    [row,Order_Max]=find(Detection==detection_Max);
    OriginalOrder_Result=OriginalOrder(Order_Max);
    [Row_Result,Col_Result]=ind2sub([num_YCoord,num_XCoord],OriginalOrder_Result);
%     row_Result=round(mean(Row_Result));
%     col_Result=round(mean(Col_Result));
%     order_Max=sub2ind([num_YCoord,num_XCoord],row_Result,col_Result);
%     if order_Max>length(Index_Min(1,:))
%         
%     end
    if length(Row_Result)==1
        Result(i_Result,:)=[Area(1,1)-0.5*gridWidth+Col_Result'*gridWidth,Area(2,2)+0.5*gridWidth-Row_Result'*gridWidth];
    else
        Result(i_Result,:)=mean([Area(1,1)-0.5*gridWidth+Col_Result'*gridWidth,Area(2,2)+0.5*gridWidth-Row_Result'*gridWidth]);
    end
    
    %擦除能量
    Origin_Index_Erasure=zeros(num_Channel,1);
    for i_Channel=1:num_Channel
        Index_MeasureOfMaxGrid=Index_Min(i_Channel,order_Max):num_Channel:Index_Max(i_Channel,order_Max);%计算最大栅格在measure中的索引
        [detection_Max,temp_Order_SingleChannelMaxEnergy]=max(Measure(Index_MeasureOfMaxGrid));%计算最大能量在Index_MeasureOfMaxGrid中对应的索引
        order_SingleChannelMaxEnergy=Index_MeasureOfMaxGrid(temp_Order_SingleChannelMaxEnergy);%计算最大能量对应在Measure中的索引
        Origin_Index_Erasure(i_Channel,1)=order_SingleChannelMaxEnergy;
    end
    Index_Erasure=kron(Origin_Index_Erasure,ones(2*num_Erasure+1,1))+repmat((-num_Erasure:num_Erasure)'*num_Channel,num_Channel,1);%计算擦除索引，按列排列
    Measure(Index_Erasure)=0;
    %在Order_Result和Index里面擦除最大的点
    Order_Result(1,order_Max)=0;
    %删除没有过门限的栅格的索引
    Index_Min=Index_Min(:,Order_Result);
    Index_Max=Index_Max(:,Order_Result);
    %删除没有过门限的栅格的原始编号
    OriginalOrder=OriginalOrder(1,Order_Result);
    
    %提取检测统计量
    num_Grid=sum(Order_Result);
    Detection=zeros(num_Channel,num_Grid);
    for i_Grid=1:num_Grid
        for i_Channel=1:num_Channel
            Detection(i_Channel,i_Grid)=max(Measure(Index_Min(i_Channel,i_Grid):num_Channel:Index_Max(i_Channel,i_Grid)));
        end
    end
    Detection=sum(Detection);
end
Result=Result(1:i_Result,:);
%------------------------------------------------------------------------
% hold on
if i_Montecarlo==1&&flag_FigureDetect==1
    plot(Result(:,1),Result(:,2),'yo','MarkerSize',10,'LineWidth',2)
    hold on
%     plot(50160,230,'yo','MarkerSize',10,'LineWidth',2)
%     hold on
    plot(Coord_Target(:,1),Coord_Target(:,2),'g+','MarkerSize',10,'LineWidth',1.5)
    title('极值搜索检测统计量平面与检测结果');
    legend('检测点迹','真实目标')
    xlabel('X/m');
    ylabel('Y/m');
    grid on
    if flag_FigureDetect_Gif
        print(gcf,'-dbmp',strcat('../Figure/Detect/Detect',num2str(i_Frame),'.bmp'));
    end
    close(temp_Figure_Detect)
end
%------------------------------------------------------------------------


