%**************************************************************************
%*��Creat time����2017-01-29 21:56          ��Version����0.0
%*��Writer����LiShuai 461837929@qq.com
%*��Writer department����
%*��Function����
%*��⣬ȥ��Ӱ
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
function Result=DetectionOfDataAlignment(Index,Measure,num_XCoord,num_YCoord,gate_Alignment,Area,Temp_XCoord,Temp_YCoord,gridWidth,num_Channel,i_Frame,Coord_Target,flag_FigureDetect,flag_FigureDetect_Gif,i_Montecarlo)
num_Erasure=2;%ǰ�������2����Ԫ
Detection=sum(Measure(Index));
OriginalOrder=(1:num_XCoord*num_YCoord);%������������ÿ��դ���ԭʼ���
Order_Result=Detection>gate_Alignment;%ÿ��������դ����Detection�еı��
num_Result=sum(Order_Result);
Result=zeros(num_Result,2);
%--------------------------------------------------------------------------
%ȥ��Ӱǰ�ĵ�
Detection_Figure=reshape(Detection,num_YCoord,num_XCoord);
[Row,Col]=find(Detection_Figure>gate_Alignment);
Coord_Result_X=Area(1,1)-0.5*gridWidth+Col*gridWidth;
Coord_Result_Y=Area(2,2)+0.5*gridWidth-Row*gridWidth;
% Temp_YCoord_Figure=zeros(size(Temp_YCoord));
% Detection_Figure=zeros(size(Detection));
% for i_Row=1:length(Temp_YCoord_Figure(:,1))
%     Temp_YCoord_Figure(i_Row,:)=Temp_YCoord(end-i_Row+1,:);
%     Detection_Figure(i_Row,:)=Detection(end-i_Row+1,:);
% end
%--------------------------------------------------------------------------
%��ͼ
if i_Montecarlo==1&&flag_FigureDetect==1
    temp_Figure_Detect=figure();
    imagesc(Temp_XCoord+gridWidth/2,Temp_YCoord-gridWidth/2,Detection_Figure)
    hold on
end
%--------------------------------------------------------------------------
% plot(50000,200,'y*','MarkerSize',10,'LineWidth',2)
% plot(Coord_Result_X,Coord_Result_Y,'g*','MarkerSize',10,'LineWidth',0.5)
%--------------------------------------------------------------------------
%ȥ��Ӱ

for i_Result=1:num_Result
    %���
    Order_Result=Detection>gate_Alignment;
    if sum(Order_Result)==0
        i_Result=i_Result-1;
        break
    end

    %�ҵ�����
    [detection_Max,order_Max]=max(Detection);
    %----------------------------------------------------------------------
    %��������
%     originalOrder_Result=OriginalOrder(order_Max);
%     [row_Result,col_Result]=ind2sub([num_YCoord,num_XCoord],originalOrder_Result);
%     Result(i_Result,:)=[Area(1,1)-0.5*gridWidth+col_Result*gridWidth,Area(2,2)+0.5*gridWidth-row_Result*gridWidth];
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %������
    [row,Order_Max]=find(Detection==detection_Max);
    OriginalOrder_Result=OriginalOrder(Order_Max);
    [Row_Result,Col_Result]=ind2sub([num_YCoord,num_XCoord],OriginalOrder_Result);
    if length(Row_Result)==1
        Result(i_Result,:)=[Area(1,1)-0.5*gridWidth+Col_Result'*gridWidth,Area(2,2)+0.5*gridWidth-Row_Result'*gridWidth];
    else
        Result(i_Result,:)=mean([Area(1,1)-0.5*gridWidth+Col_Result'*gridWidth,Area(2,2)+0.5*gridWidth-Row_Result'*gridWidth]);
    end
    %----------------------------------------------------------------------
    %��������
    Origin_Index_Erasure=Index(:,order_Max);%��������������
    Index_Erasure=kron(Origin_Index_Erasure,ones(2*num_Erasure+1,1))+repmat((-num_Erasure:num_Erasure)'*num_Channel,num_Channel,1);%���������������������
    Measure(Index_Erasure)=0;
    %��Order_Result��Index����������ĵ�
    Order_Result(1,order_Max)=0;
    %ɾ��û�й����޵�դ�������
    Index=Index(:,Order_Result);
    %ɾ��û�й����޵�դ���ԭʼ���
    OriginalOrder=OriginalOrder(1,Order_Result);
    Detection=sum(Measure(Index));
end
Result=Result(1:i_Result,:);
%--------------------------------------------------------------------------
% ��ͼ
if i_Montecarlo==1&&flag_FigureDetect==1
    plot(Result(:,1),Result(:,2),'yo','MarkerSize',10,'LineWidth',2)
    hold on
    plot(Coord_Target(:,1),Coord_Target(:,2),'g+','MarkerSize',10,'LineWidth',1.5)
    title('���ݶ�����ͳ����ƽ��������');
    xlabel('X/m');
    ylabel('Y/m');
    legend('���㼣','��ʵĿ��')
    grid on
    if flag_FigureDetect_Gif
        print(gcf,'-dbmp',strcat('../Figure/Detect/Detect',num2str(i_Frame),'.bmp'));
    end
    close(temp_Figure_Detect)
end
%--------------------------------------------------------------------------
%%
%��ͼ
% print(gcf,'-dbmp',strcat('../Figure/Detection',num2str(i_Frame),'.bmp'));


