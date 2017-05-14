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

%% �˲������£�KF��  
 Xk = F1*X; % һ��Ԥ��״̬����
 Zk = H*Xk; % Ԥ������ֵ���˲����������ģ�
 Pk = F1*P*F1' + Q1; % һ��Ԥ�����Э�������
 S = H*Pk*H' + R; % ��Ϣ��������ؾ���
 K = (Pk*H')*inv(S); % �˲�������

 %% ѡȡ��Ч����ֵ
 switch process,
    case 'Non_AI',        
        vv  = zeros(2,length(x)); dv  = zeros(1,length(x));        
        V = zeros(2,1); % ��Ϣ 
        for ii = 1:length(x)
             vv(:,ii)  = [x(ii);y(ii)]-Zk; 
             dv(ii)  = vv(:,ii)'*inv(S)*vv(:,ii); 
             if dv(ii) > gama
                 dv(ii)  = 0;          
             end
        end 
        index = find(dv); mk = length(index); % �����ڵ�������
        if any(mk) == 1
            vv_in = vv(:,index); % ���벨���ڵ���Ч�������Ϣֵ
            dv_in = dv(index); % ���ӿ����ֲ�    
            beta = zeros(1,length(index)+1); % ��������
            e = exp(-0.5*dv_in);      
            b = mk*(1-Pd*Pg)/(Pd*2*gama*(det(S))^0.5); % ���� b �ǲ�����PDA
            beta(1) = b/(b+sum(e));  beta(2:end) = e./(b+sum(e)); % beta0��beta1-mk      
            V(1) = beta(2:end)*vv_in(1,:)';
            V(2) = beta(2:end)*vv_in(2,:)'; % ������Ϣ���̣�V(k)=sum(beta1(i)*V(i)(k))��
        end
     case 'AI',
         vv  = zeros(2,length(x)); dv  = zeros(1,length(x));        
         V = zeros(2,1); % ��Ϣ 
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
         index = find(dv); mk = length(index); % �����ڵ�������
         if any(mk) == 1
             vv_in = vv(:,index); % ���벨���ڵ���Ч�������Ϣֵ
             amp_r = ratio(index);
             dv_in = dv(index); % ���ӿ����ֲ�    
             beta = zeros(1,length(index)+1); % ��������
             e = exp(-0.5*dv_in);      
             b = mk*(1-Pd*Pg)/(Pd*2*gama*(det(S))^0.5); % ���� b �ǲ�����PDA
             beta(1) = b/(b+sum(e*amp_r'));  beta(2:end) = (e.amp_r)./(b+sum(e*amp_r')); % beta0��beta1-mk      
             V(1) = beta(2:end)*vv_in(1,:)';
             V(2) = beta(2:end)*vv_in(2,:)'; % ������Ϣ���̣�V(k)=sum(beta1(i)*V(i)(k))��
         end
 end

%% �������˲�����״̬���� 
if any(mk) == 1
    X = Xk+K*V; % ״̬���£�x(k|k)=x(k|k-1)+W(k)*V(k)��
    flag_Inside=1;
else
    X = Xk; 
end

%% �������Э�������
if any(mk) == 1
    vw=zeros(2);
    for ii = 1:mk
         vw=vw+beta(ii+1)*vv_in(:,ii)*vv_in(:,ii)';
    end
    P_hat = K*(vw-V*V')*K';
    Pc = (eye(4)-K*H)*Pk;
    P = Pk*beta(1)+(1-beta(1))*Pc+P_hat;
else
    ct = (1-9*exp(-8))/(1-exp(-8)); % ����b0ʱ��һ�����������걸٤����
    b0 = Pd*Pg*(1-ct)/(1-Pd*Pg); % Сֵ����
    P = Pk+b0*K*S*K';
end

%% ���
Xpda = X; P_estimate = P;
