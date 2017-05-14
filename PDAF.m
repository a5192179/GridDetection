function [Xpda,P_estimate,flag_Inside] = PDAF(Result_Detection,Xpda,P_estimate,mode)
flag_Inside=0;
%% Input Judgment
if nargin<3 || nargin>5
    error( 'Wrong number of input arguments');
end
if mode == 1
    process = 'AI';
elseif mode == 0
    process = 'Non_AI';
end

%% Input Parameters Distribution
global F1 H Q1 R SNR Pfa Pd Pg gama
% [x,y] = find(data_detection); % the location of radar return on the Data Plane per frame
x=Result_Detection(:,1);
y=Result_Detection(:,2);
P = P_estimate; X = Xpda;

%% 滤波器更新（KF）  
 Xk = F1*X; % 一步预测状态向量
 Zk = H*Xk; % 预测量测值（滤波器波门中心）
 Pk = F1*P*F1' + Q1; % 一步预测误差协方差矩阵
 S = H*Pk*H' + R; % 新息过程自相关矩阵
 K = (Pk*H')*inv(S); % 滤波器增益

 %% 选取有效量测值
 switch process,
    case 'Non_AI',        
        vv  = zeros(2,length(x)); dv  = zeros(1,length(x));        
        V = zeros(2,1); % 新息 
        for ii = 1:length(x)
             vv(:,ii)  = [x(ii);y(ii)]-Zk; 
             dv(ii)  = vv(:,ii)'*inv(S)*vv(:,ii); 
             if dv(ii) > gama
                 dv(ii)  = 0;          
             end
        end 
        index = find(dv); mk = length(index); % 波门内的量测数
        if any(mk) == 1
            vv_in = vv(:,index); % 落入波门内的有效量测的信息值
            dv_in = dv(index); % 服从卡方分布    
            beta = zeros(1,length(index)+1); % 关联概率
            e = exp(-0.5*dv_in);      
            b = mk*(1-Pd*Pg)/(Pd*2*gama*(det(S))^0.5); % 参数 b 非参数化PDA
            beta(1) = b/(b+sum(e));  beta(2:end) = e./(b+sum(e)); % beta0、beta1-mk      
            V(1) = beta(2:end)*vv_in(1,:)';
            V(2) = beta(2:end)*vv_in(2,:)'; % 联合新息过程（V(k)=sum(beta1(i)*V(i)(k))）
        end
     case 'AI',
         vv  = zeros(2,length(x)); dv  = zeros(1,length(x));        
         V = zeros(2,1); % 新息 
         ratio = zeros(1,length(x));
         for ii = 1:length(x)
              vv(:,ii)  = [x(ii);y(ii)]-Zk; 
              dv(ii)  = vv(:,ii)'*inv(S)*vv(:,ii);
              ratio(ii) = Pfa/(Pd*(1+10^(SNR/10)))*exp(data_detection(x(ii),y(ii))^2/(2*(1-1/(1+10^(SNR/10)))));
              if dv(ii) > gama
                  dv(ii)  = 0;
                  ratio(ii) = 1;
              end
         end 
         index = find(dv); mk = length(index); % 波门内的量测数
         if any(mk) == 1
             vv_in = vv(:,index); % 落入波门内的有效量测的信息值
             amp_r = ratio(index);
             dv_in = dv(index); % 服从卡方分布    
             beta = zeros(1,length(index)+1); % 关联概率
             e = exp(-0.5*dv_in);      
             b = mk*(1-Pd*Pg)/(Pd*2*gama*(det(S))^0.5); % 参数 b 非参数化PDA
             beta(1) = b/(b+sum(e*amp_r'));  beta(2:end) = (e.amp_r)./(b+sum(e*amp_r')); % beta0、beta1-mk      
             V(1) = beta(2:end)*vv_in(1,:)';
             V(2) = beta(2:end)*vv_in(2,:)'; % 联合新息过程（V(k)=sum(beta1(i)*V(i)(k))）
         end
 end

%% 更新子滤波器的状态变量 
if any(mk) == 1
    X = Xk+K*V; % 状态更新（x(k|k)=x(k|k-1)+W(k)*V(k)）
    flag_Inside=1;
else
    X = Xk; 
end

%% 更新误差协方差矩阵
if any(mk) == 1
    vw=zeros(2);
    for ii = 1:mk
         vw=vw+beta(ii+1)*vv_in(:,ii)*vv_in(:,ii)';
    end
    P_hat = K*(vw-V*V')*K';
    Pc = (eye(4)-K*H)*Pk;
    P = Pk*beta(1)+(1-beta(1))*Pc+P_hat;
else
    ct = (1-9*exp(-8))/(1-exp(-8)); % 计算b0时的一个参数，不完备伽马函数
    b0 = Pd*Pg*(1-ct)/(1-Pd*Pg); % 小值参数
    P = Pk+b0*K*S*K';
end

%% 输出
Xpda = X; P_estimate = P;
