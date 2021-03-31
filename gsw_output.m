%% gsw_output -  Calculate the phase pattern to generate multiple spot pattern using phase-only SLM.
% The input parameters: size_SLM is the size of SLM, which is presented as
% a vector with two elements. e.g. size_SLM=[1920,1080] for a 1k resolution
% weight is the beam array intensity, for example one(10,10) for a 10x10
% beam array with uniform intensity
% interval is the distance between the arrays compared to the beam size
% [Image_SLM,Error1,Error2]=gsw_output(size_SLM,weight,interval) returns
% the phase pattern encoded on SLM and the RMSE of error between target
% weight and simulated weight during each iterations.


function [phase_SLM,Error1,Error2] = gsw_output(size_SLM,weight,interval)
%% This function is defined to calculate two type error:
%if type==0, then the whole area is calculated. if type==1, then using
%small area to calculate
    function [error,Ik]=ErrorCal(phi,position,type)
        size_=(size_FFT-1)/2;
        [X,Y] = meshgrid(-size_(1):1:size_(1), -size_(2):1:size_(2));
        A0=exp( - ((X').^2)/(1000^2) - ((Y').^2)/(1000^2) ).*exp(1i*phi);
        if type==1
            A0=A0(position_time(1,1):position_time(1,2),position_time(2,1):position_time(2,2));
        end
        B=fftshift(fft2(A0,size_FFT(1),size_FFT(2)));
        Ik=IntensityMeasure(B,position);
        error=(Ik/mean(Ik(:)))-(weight/mean(weight(:)));
    end

%% This function is defined to measure the simulation beam power
    function [power,power_sum]=IntensityMeasure(B,position)
        rc=(position(:,2)-position(:,1)+1)/interval;
        power=zeros(rc(1),rc(2));
        for i=1:rc(1)
            for ii=1:rc(2)
                x=position(1,1)+interval/2+(i-1)*interval;
                y=position(2,1)+interval/2+(ii-1)*interval;
                ratio_x=ceil(ratio*size_SLM(1)/size_SLM(2));
                power(i,ii)=sum(sum(abs(B(x-ratio_x/2:x+ratio_x/2-1,y-ratio/2:y+ratio/2-1)).^2));
            end
        end
        power_sum=sum(abs(B(:)).^2);
    end
%% This function is the algorithm for GS
    function [Output,g_next]=GS_algorithm(phase,g,position)
    %Step1: Gaussian beam with given phase    
        size_=(size_FFT-1)/2;
        [X,Y] = meshgrid(-size_(1):1:size_(1), -size_(2):1:size_(2));
        A0=exp( - ((X').^2)/(1000^2) - ((Y').^2)/(1000^2) ).*exp(1i*phase);

        B0=fftshift(fft2(A0,size_FFT(1),size_FFT(2)));
        A0=A0(position_time(1,1):position_time(1,2),position_time(2,1):position_time(2,2));
        %Step2: fft the beam
        
        B=fftshift(fft2(A0,size_FFT(1),size_FFT(2)));
        
        ak=sqrt(IntensityMeasure(B,position));
        
        g_next=(sqrt(weight)/sum(sqrt(weight(:))))./(ak/sum(ak(:))).*g;
        
        [at,~]=Multibeam(g_next.*sqrt(weight));
        
        %Step3: Replace the amplitude by given intensity
        D=(at).*exp(1i*angle(B0));
        %Step4: ifft get new phase
        E=ifft2(ifftshift(D));
        Output=angle(E);
        
%         phase_small=zeros(size(phase));
%         phase_small(position(1,1):position(1,2),position(2,1):position(2,2))=phase(position(1,1):position(1,2),position(2,1):position(2,2));
%         A0_small=exp( - ((X').^2)/(500^2) - ((Y').^2)/(500^2) ).*exp(1i*phase_small);
%         B_small=fftshift(fft2(A0_small,size_part(1),size_part(2)));
%         power=sqrt(IntensityMeasure(B_small,position))
%         power.^2/sum(power(:).^2)-weight/sum(weight(:))
    end
%% This function is designed for generate pattern of arbitary multibeam
    function [Multipattern,position]=Multibeam(weight)
        [row, column]=size(weight);
        
        single_r=(interval-1)/2;
        [single_x,single_y]=meshgrid([-single_r:single_r],[-single_r:single_r]);
        singlepattern=exp(-2*(single_x.^2+single_y.^2)/w0^2);
        Multi=repmat(singlepattern,row,column);

        for i=1:row
            for ii=1:column
                Multi((i-1)*interval+1:i*interval,(ii-1)*interval+1:(ii)*interval)=Multi((i-1)*interval+1:i*interval,(ii-1)*interval+1:(ii)*interval)*weight(i,ii);
            end
        end     
        
        [Multi_x,Multi_y]=size(Multi);
        Multipattern=zeros(size_FFT);       
        position=[floor(size_FFT(1)/2)-floor(Multi_x/2)+1,floor(size_FFT(1)/2)+floor(Multi_x/2);floor(size_FFT(2)/2)-floor(Multi_y/2)+1,floor(size_FFT(2)/2)+floor(Multi_y/2)];
        Multipattern(position(1,1):position(1,2),position(2,1):position(2,2))=Multi;   
    end
% %% This function can generate grating 
     function grate=grating(size,grating_Tx,grating_Ty)
        grate_x=mod(repmat([1:size(1)]',1,size(2)),grating_Tx)/grating_Tx;
        grate_y=mod(repmat([1:size(2)],size(1),1),grating_Ty)/grating_Ty;
        grate=mod(grate_x+grate_y,1)*2*pi;
     end
% %% This function can generate modulation
%     function phase_mod=phase_mod(size,mod_T,mod_depth)
%         num=floor(size(1)/abs(mod_T));
%         phase_mod= mod_depth*pi*repmat(sin([1:mod_T]'/mod_T*2*pi),num,size(2));
%     end

%% Main function
tic

%parameter for target beam, w0 is the waist of reflected beams. 'ratio'
%here is to enlarge the size of time domin of FFT to achieve high precision in
%frequency domain. Noticed that larger the ratio, output more precise but more time needed
w0=1;
ratio=2;

size_FFT=[1 1]*size_SLM(1)*ratio;
padnum=(size_FFT-size_SLM)./2;
position_time=[padnum(1)+1,padnum(1)+size_SLM(1);padnum(2)+1,padnum(2)+size_SLM(2)];

[~,position_freq]=Multibeam(sqrt(weight)); %intensity

%Randomly generate initial phase
Phase0=rand(size_FFT);
g=ones(size(weight));

[phi,~]=GS_algorithm(Phase0,g,position_freq);
[Error1(:,:,1),Ik1(:,:,1)]=ErrorCal(phi,position_freq,0);
[Error2(:,:,1),Ik2(:,:,1)]=ErrorCal(phi,position_freq,1);

% The iteration time is set to be 10, which can be adjust by user
for nn=1:10
    nn
    [phi,g]=GS_algorithm(phi,g,position_freq);
    [Error1(:,:,nn+1),Ik1(:,:,nn+1)]=ErrorCal(phi,position_freq,0);
    [Error2(:,:,nn+1),Ik2(:,:,nn+1)]=ErrorCal(phi,position_freq,1);
end

Phase_f=phi(position_time(1,1):position_time(1,2),position_time(2,1):position_time(2,2));

% User can also add grating like below to adjust the position of beam array
%Phase_f=Phase_f+grating(fliplr(size_real),10)'+grating(size_real,-10);

Eout=exp(1i*Phase_f);   
Eout=fftshift(fft2(Eout,size_FFT(1),size_FFT(2)));
Iout=abs(Eout).^2;
imshow(I_out/max(Iout))

phase_SLM=mod(Phase_all,2*pi);
% Image_SLM for user to export Image uint8 format  
% Image_SLM=uint8(Phase_n*255/(2*3.1416));
% imshow(Image_SLM);
toc

end