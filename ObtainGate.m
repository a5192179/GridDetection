%**************************************************************************
%*【Creat time】：2017-02-01 12:43          【Version】：0.0
%*【Writer】：LiShuai 461837929@qq.com
%*【Writer department】：
%*【Function】：
%*获得门限
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
clc
clear all
close all
time=datestr(now);
fprintf('Begin time:%s\n',time)
%参数设置
num_channel=4;
saveRootFolder='..\Data';
saveFileFolder=[saveRootFolder,'\Gate'];
if ~exist(saveFileFolder,'dir') 
    mkdir(saveFileFolder)
end 
saveFileName_Alignment=[saveFileFolder,'\AlignmentGate',num2str(num_channel,'%02d'),'.txt'];
saveFileName_Elect=[saveFileFolder,'\ElectGate',num2str(num_channel,'%02d'),'.txt'];
pd=1e-5;
num_MonteCarlo=1e7;
%%
%对齐的门限
Energy_Alignment=sort(sum((abs(sqrt(1/2)*(randn(num_MonteCarlo,num_channel)+1i*randn(num_MonteCarlo,num_channel)))).^2,2));
gate_Alignment=Energy_Alignment(end-(num_MonteCarlo*pd-1));
save(saveFileName_Alignment,'gate_Alignment','-ascii')
fprintf('Alignment gate=%d\n',gate_Alignment)
%%
% Length_Sig=(1:30)';
% Gate_Elect=zeros(length(Length_Sig(:,1)),1);
% for i_Length_Sig=1:length(Length_Sig(:,1))
%     length_Sig=Length_Sig(i_Length_Sig,1);
%     Energy_Elect=zeros(1,num_MonteCarlo);
%     for i_channel=1:num_channel
%         Noise=sqrt(1/2)*(randn(length_Sig,num_MonteCarlo)+1i*randn(length_Sig,num_MonteCarlo));
%         if length_Sig~=1
%             Energy_Elect=Energy_Elect+max((abs(Noise)).^2);
%         else
%             Energy_Elect=Energy_Elect+(abs(Noise)).^2;
%         end
%     end
%     Energy_Elect=sort(Energy_Elect);
% %     Energy_Elect=Energy_Elect';
%     Gate_Elect(i_Length_Sig,1)=Energy_Elect(end-(num_MonteCarlo*pd-1));
%     fprintf('Elect gate=%d,time:%s\n',Gate_Elect(i_Length_Sig,1),time)
%     time=datestr(now);
% end
% save(saveFileName_Elect,'Gate_Elect','-ascii')
% fprintf('End time:%s\n',time)
