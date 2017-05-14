function [data_Alarm] = CFAR(data,N_Azi,N_APro,N_Range,N_RPro,Pfa,mode)
% This function is applied to process measured data by traditional detection method！！CFAR
% % All of the data have the same threshold, so mean(sum(cells))
% data_Alarm: The result CFAR detection. The scale of this matrix is as large as data, 
%             if some cell of data surpass the threshold the cell of data_Alarm would be one, if else whould be zero.
% data   :The data plane . It should be a one or two dimension data. Row means Azimuth, column means Range
% N_Range:The length of sliding window along Range direction (include detection cell and protected cells)
%         it always is odd and 0<N_Range<Range_max
% N_RPro :The length of protected cells along Range direction(just 1 or 0)
% N_Azi  :The length of sliding window along Azimuth direction(include detection cell and protected cells)
%         it always is odd 0<N_Azi<Azimuth_max
% N_APro :The length of protected cells along Azimuth direction
% Pfa    :The false alarm rate 
% mode   :The tow dimension CFAR mode for the measured data(rectangle<default> or cross)
% data_Alarm=CFAR_measured_data(data,N_Azi,N_APro,N_Range,N_RPro,T,'rectangle')
if nargin<6 || nargin>7
    error( 'Wrong number of input arguments');
elseif nargin==6
    mode='rectangle';
end
if ~(strcmp(mode,'rectangle') || strcmp(mode,'cross'))
    error( 'Wrong: mode should be ''rectangle'' or ''cross''.'); 
end
if mod(N_Range,2)~=1 || mod(N_Azi,2)~=1 || N_APro>1 || N_RPro>1 || N_APro<0 || N_RPro<0 
    error('Wrong: N_Range or N_Azi should be odd integer,and N_APro or N_RPro should be integer'); 
end
if N_Range<=2*N_RPro+1 && N_Azi<=2*N_APro+1, error('Wrong: N_Range<=2*N_RPro+1 and N_Azi<=2*N_APro+1'); end
data=abs(data).^2; % Amplitude
[NA,NR]=size(data);
if N_Range>NR || N_Azi>NA, error('Wrong: length of Range windle should be less than data column NR and N_Azi>NA (Azimuth_max)'); end
data_Alarm=zeros(NA,NR);
% tt = []; % 朔紗議
% 
% % The edge of data is divided from 1、2、3、4、6、7、8、9 areas, each area has unique processing method.
% 
%         |------------------->Range      
%         |-1-|-2-|-3-|
% data    |-4-|-5-|-6-|
%         |-7-|-8-|-9-|
%         |
%    Azimuth
% 
switch mode,
    case 'rectangle',
    for Ri=1:NR
        for Ai=1:NA
%           area 5
            if Ai-floor(N_Azi/2)>0 && Ai+floor(N_Azi/2)<=NA && Ri-floor(N_Range/2)>0 && Ri+floor(N_Range/2)<=NR
              data0=sum(sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri-floor(N_Range/2):Ri+floor(N_Range/2))))...
                  -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
              N_aver = (N_Azi*N_Range-(2*N_APro+1)*(2*N_RPro+1));
              data0=data0/N_aver;
              T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
%           area 1    
            elseif Ai-floor(N_Azi/2)<=0 && Ri-floor(N_Range/2)<=0
                if Ai==1 && Ri==1
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(1:Ai+N_APro,1:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*ceil(N_Range/2)-(N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==1 && Ri~=1
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(1:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*(ceil(N_Range/2)+Ri-1)-(N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=1 && Ri==1
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,1:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)*ceil(N_Range/2)-(2*N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(sum(data(1:Ai+floor(N_Azi/2),1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)*(ceil(N_Range/2)+Ri-1)-(2*N_APro+1)*(2*N_RPro+1)); 
                    data0=data0/N_aver; 
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)<=0 && Ri-floor(N_Range/2)>0 && Ri+floor(N_Range/2)<=NR
%           area 2              
                if Ai==1
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),Ri-floor(N_Range/2):Ri+floor(N_Range/2))))...
                        -sum(sum(data(1:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*N_Range-(N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),Ri-floor(N_Range/2):Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)*N_Range-(2*N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)<=0 && Ri+floor(N_Range/2)>NR
%           area 3 
                if Ai==1 && Ri==NR
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(1:Ai+N_APro,Ri-N_RPro:NR)));
                    N_aver = (ceil(N_Azi/2)*ceil(N_Range/2)-(N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==1 && Ri~=NR
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(1:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*(ceil(N_Range/2)+NR-Ri)-(N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=1 && Ri==NR
                    data0=sum(sum(data(1:Ai+floor(N_Azi/2),Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:NR)));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)*ceil(N_Range/2)-(2*N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(sum(data(1:Ai+floor(N_Azi/2),Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)*(ceil(N_Range/2)+NR-Ri)-(2*N_APro+1)*(2*N_RPro+1)); 
                   data0=data0/N_aver;
                   T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)>0 && Ai+floor(N_Azi/2)<=NA && Ri-floor(N_Range/2)<=0
%           area 4
                if Ri==1
                    data0=sum(sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,1:Ri+N_RPro)));
                    N_aver = (N_Azi*ceil(N_Range/2)-(2*N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (N_Azi*(ceil(N_Range/2)+Ri-1)-(2*N_APro+1)*(2*N_RPro+1)); 
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)>0 && Ai+floor(N_Azi/2)<=NA && Ri+floor(N_Range/2)>NR
%           area 6
                if Ri==NR
                    data0=sum(sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:NR)));
                    N_aver = (N_Azi*ceil(N_Range/2)-(2*N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (N_Azi*(ceil(N_Range/2)+NR-Ri)-(2*N_APro+1)*(2*N_RPro+1)); 
                    data0=data0/N_aver; 
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
%           area 7    
            elseif Ai+floor(N_Azi/2)>NA && Ri-floor(N_Range/2)<=0
                if Ai==NA && Ri==1
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:NA,1:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*ceil(N_Range/2)-(N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==NA && Ri~=1
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:NA,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*(ceil(N_Range/2)+Ri-1)-(N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=NA && Ri==1
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,1:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)*ceil(N_Range/2)-(2*N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(sum(data(Ai-floor(N_Azi/2):NA,1:Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)*(ceil(N_Range/2)+Ri-1)-(2*N_APro+1)*(2*N_RPro+1)); 
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end                
            elseif Ai+floor(N_Azi/2)>NA && Ri-floor(N_Range/2)>0 && Ri+floor(N_Range/2)<=NR
%           area 8              
                if Ai==NA
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,Ri-floor(N_Range/2):Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:NA,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*N_Range-(N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,Ri-floor(N_Range/2):Ri+floor(N_Range/2))))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)*N_Range-(2*N_APro+1)*(2*N_RPro+1)); 
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai+floor(N_Azi/2)>NA && Ri+floor(N_Range/2)>NR
%           area 9 
                if Ai==NA && Ri==NR
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:NA,Ri-N_RPro:NR)));
                    N_aver = (ceil(N_Azi/2)*ceil(N_Range/2)-(N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==NA && Ri~=NR
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:NA,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = (ceil(N_Azi/2)*(ceil(N_Range/2)+NR-Ri)-(N_APro+1)*(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=NA && Ri==NR
                    data0=sum(sum(data(Ai-floor(N_Azi/2):NA,Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:NR)));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)*ceil(N_Range/2)-(2*N_APro+1)*(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(sum(data(Ai-floor(N_Azi/2):NA,Ri-floor(N_Range/2):NR)))...
                        -sum(sum(data(Ai-N_APro:Ai+N_APro,Ri-N_RPro:Ri+N_RPro)));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)*(ceil(N_Range/2)+NR-Ri)-(2*N_APro+1)*(2*N_RPro+1));
                   data0=data0/N_aver; 
                   T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            else
                
            end
%             tt = [tt sqrt(T*data0)]; % check the threashold
    if data(Ai,Ri)>T*data0
%         data_Alarm(Ai,Ri)=1;
          data_Alarm(Ai,Ri) = sqrt(data(Ai,Ri));
    end
        end
    end
% one of two dimension CA-CFAR methods！！cross data
    case 'cross',
    for Ri=1:NR
        for Ai=1:NA
%           area 5
            if Ai-floor(N_Azi/2)>0 && Ai+floor(N_Azi/2)<=NA && Ri-floor(N_Range/2)>0 && Ri+floor(N_Range/2)<=NR
              data0=sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri))...
                  +sum(data(Ai,Ri-floor(N_Range/2):Ri+floor(N_Range/2)))...
                  -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
              N_aver = (N_Azi+N_Range-(2*N_APro+1)-(2*N_RPro+1));
              data0=data0/N_aver;
              T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
%           area 1    
            elseif Ai-floor(N_Azi/2)<=0 && Ri-floor(N_Range/2)<=0
                if Ai==1 && Ri==1
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(1:Ai+N_APro,Ri))-sum(data(Ai,1:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+ceil(N_Range/2)-(N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==1 && Ri~=1
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(1:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+(ceil(N_Range/2)+Ri-1)-(N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=1 && Ri==1
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,1:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)+ceil(N_Range/2)-(2*N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                       +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)+(ceil(N_Range/2)+Ri-1)-(2*N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)<=0 && Ri-floor(N_Range/2)>0 && Ri+floor(N_Range/2)<=NR
%           area 2              
                if Ai==1
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):Ri+floor(N_Range/2)))...
                        -sum(data(1:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+N_Range-(N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)+N_Range-(2*N_APro+1)-(2*N_RPro+1)); 
                    data0=data0/N_aver; 
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)<=0 && Ri+floor(N_Range/2)>NR
%           area 3 
                if Ai==1 && Ri==NR
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(1:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:NR));
                    N_aver = (ceil(N_Azi/2)+ceil(N_Range/2)-(N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==1 && Ri~=NR
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(1:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+(ceil(N_Range/2)+NR-Ri)-(N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=1 && Ri==NR
                    data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:NR));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)+ceil(N_Range/2)-(2*N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(data(1:Ai+floor(N_Azi/2),Ri))...
                       +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+Ai-1)+(ceil(N_Range/2)+NR-Ri)-(2*N_APro+1)-(2*N_RPro+1)); 
                   data0=data0/N_aver; 
                   T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)>0 && Ai+floor(N_Azi/2)<=NA && Ri-floor(N_Range/2)<=0
%           area 4
                if Ri==1
                    data0=sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,1:Ri+N_RPro));
                    N_aver = (N_Azi+ceil(N_Range/2)-(2*N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (N_Azi+(ceil(N_Range/2)+Ri-1)-(2*N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai-floor(N_Azi/2)>0 && Ai+floor(N_Azi/2)<=NA && Ri+floor(N_Range/2)>NR
%           area 6
                if Ri==NR
                    data0=sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:NR));
                    N_aver = (N_Azi+ceil(N_Range/2)-(2*N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(data(Ai-floor(N_Azi/2):Ai+floor(N_Azi/2),Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (N_Azi+(ceil(N_Range/2)+NR-Ri)-(2*N_APro+1)-(2*N_RPro+1)); 
                    data0=data0/N_aver; 
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
%           area 7    
            elseif Ai+floor(N_Azi/2)>NA && Ri-floor(N_Range/2)<=0
                if Ai==NA && Ri==1
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:NA,Ri))-sum(data(Ai,1:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+ceil(N_Range/2)-(N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==NA && Ri~=1
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:NA,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+(ceil(N_Range/2)+Ri-1)-(N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=NA && Ri==1
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,1:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)+ceil(N_Range/2)-(2*N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                       +sum(data(Ai,1:Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)+(ceil(N_Range/2)+Ri-1)-(2*N_APro+1)-(2*N_RPro+1)); 
                    data0=data0/N_aver; 
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end                
            elseif Ai+floor(N_Azi/2)>NA && Ri-floor(N_Range/2)>0 && Ri+floor(N_Range/2)<=NR
%           area 8              
                if Ai==NA
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:NA,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+N_Range-(N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):Ri+floor(N_Range/2)))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)+N_Range-(2*N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            elseif Ai+floor(N_Azi/2)>NA && Ri+floor(N_Range/2)>NR
%           area 9 
                if Ai==NA && Ri==NR
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:NA,Ri))-sum(data(Ai,Ri-N_RPro:NR));
                    N_aver = (ceil(N_Azi/2)+ceil(N_Range/2)-(N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai==NA && Ri~=NR
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:NA,Ri))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = (ceil(N_Azi/2)+(ceil(N_Range/2)+NR-Ri)-(N_APro+1)-(2*N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                elseif Ai~=NA && Ri==NR
                    data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                        +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:Ai+N_APro,Ri))-sum(data(Ai,Ri-N_RPro:NR));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)+ceil(N_Range/2)-(2*N_APro+1)-(N_RPro+1));
                    data0=data0/N_aver;
                    T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                else
                   data0=sum(data(Ai-floor(N_Azi/2):NA,Ri))...
                       +sum(data(Ai,Ri-floor(N_Range/2):NR))...
                        -sum(data(Ai-N_APro:Ai+N_APro))-sum(data(Ai,Ri-N_RPro:Ri+N_RPro));
                    N_aver = ((ceil(N_Azi/2)+NA-Ai)+(ceil(N_Range/2)+NR-Ri)-(2*N_APro+1)-(2*N_RPro+1)); 
                   data0=data0/N_aver; 
                   T = N_aver*(Pfa^(-1/N_aver) - 1); % multiplier factor
                end
            else
                
            end
    if data(Ai,Ri)>T*data0
%         data_Alarm(Ai,Ri)=1;
          data_Alarm(Ai,Ri) = sqrt(data(Ai,Ri));
    end
        end
    end
end