clear all
clc
%Load data
% data_1= struct2array(load('data1_rec1208-145650_r@441Hz.mat')); 
ENF= struct2array(load('ENFdata_FDR.mat'));
data_1 = struct2array(load('rec1208-145644_r@441Hz.mat')); %% this is data 2
%% signal filtering 
Fs= 441;                % sampling frequency (Hz)
Ts= 1/Fs;               % sampling period (seconds)
F_Nyq=Fs/2;             % Nyquist frequency (Hz)
F_pass1=119.95;          %First cut-off frequency (Hz)
F_pass2=120.05;          %Second cut-off frequency (Hz)
F_pass1_n=F_pass1/F_Nyq;     %Normalized first cut-off frequency
F_pass2_n=F_pass2/F_Nyq;     %Normalized second cut-off frequency

N= 1501;                %Filter order (must be odd number)
win=hamming(N);         %Window function

b= fir1(N-1,[F_pass1_n F_pass2_n],'bandpass',win,'scale'); 
   %Filter coefficients (impulse response)
data_filtered1 = filtfilt(b,1,data_1); %Apply filtering

clear b
clear F_Nyq
clear F_pass1
clear F_pass1_n
clear F_pass2
clear F_pass2_n
clear N
clear win
clear Ts
%% Spectral estimation
Fs = 441;
frame_len_secs_data1=5; 
frame_length= frame_len_secs_data1 * Fs; % frame length in samples
M=frame_length;
nfft= 2*M;
K=nfft;
signal_len=length(data_1);          
shift_amount= Fs;  
w1 = rectwin(M);
rown = ceil((1+K)/2); % calculate the total number of rows
coln = 1+fix((signal_len-M)/shift_amount);    
threshold=1e-6;

k=0:K-1; 
M_vec=0:M-1;
omega_f=2*pi*k/K;
sigma2=1e-8;

break_time=0;
indx = 0;
col = 1;
tic
while indx + M <= (signal_len/2)+Fs;
    value=20e10;
     disp(['frame segment: ', num2str(col)]);
    
    % initialization
    xw = data_filtered1(indx+1:indx+M).*w1;
    iaa_initial_est = ifft(xw,K)*(K/M);
    p=abs(iaa_initial_est).^2; %K X 1
    newp=p;
    while value>threshold
       rr=fft(p,K);
       %locate the elements of the first column of R
       r_M=rr(1:M,1);
       r_M(1,1)=r_M(1,1)+sigma2; %for numerical stability
       [m1,m2,m3,m4]=topinv(conj(r_M),conj(r_M)');
       %Generators of Rinv
       t_M=m1{1,1}; %column vector
       w_Mminus1=t_M(2:end);
       mm2=m2{1,2}; 
       %error2=1/mm2(1,1); %first element of the row vector=1/sigma^2
       t_M=t_M*sqrt(mm2(1,1));
       %s_M=m3{1,1}; %column vector
       %s_M=s_M*sqrt(mm2(1,1));
       s_M=vertcat(0,conj(flipud(w_Mminus1)))*sqrt(mm2(1,1));
       
       %debug from this point onwards
       I_Mminus1=eye(M-1);
       Z_M=zeros(M,M);
       Z_M(2:end,1:end-1) = I_Mminus1; % unit lower shift matrix 
       Z_M2 = Z_M;
       Z_M2(1,M) = 1;
       I=eye(M);
       J = fliplr(I); % exchange matrix
       
       gg1=gallery('krylov',Z_M,t_M,M);
       gg2=gallery('krylov',Z_M,s_M,M);
       R_inv=gg1*gg1'-gg2*gg2';
       
       %Computation of the polynomial in the numerator of IAA
       zz_M=R_inv*xw;
       numerator = ifft(zz_M,K)*K;
       % Computation of the polynomial in the denominator of IAA
       vv=1:1:M;
       tt=flipud(t_M);
       tt_M=vv'.* tt; %\tilde{t_M}
       ss=flipud(s_M);
       ss_M=vv'.* ss; %\tilde{s_M}
       cc1=gallery('krylov',Z_M,tt_M,M)*conj(t_M);
       cc2=gallery('krylov',Z_M,ss_M,M)*conj(s_M);
       cc_=cc1-cc2; %\underline{c}
       cc=vertcat(cc_(M), flipud(conj(cc_(1:M-1))),zeros(K-2*M+1,1),cc_(1:M-1));
       denominator=fft(cc,K);
   
       newp_temp = numerator./denominator;
       newp=abs(newp_temp).^2;
                   
       value=((newp-p)'*(newp-p))/(p'*p);
        disp(['criterion: ', num2str(value)]);
       p=newp; 
    end                                       
    iaa_sq=p;
    
    %r(:,col) = snr(p);    
    %update the TRIAA matrix
    triaa(:,col)=iaa_sq(1:rown);
   % save temp_triaa.mat triaa
    % update the indexes
    indx = indx + shift_amount;
    %tic
    col = col + 1;
    %toc 
end
toc

triaa = triaa/M;

%triaa = triaa;
power_vector=log10((triaa));  % compute power vector
[~, index] = max(power_vector);       % find indices of max power
f = (0:rown-1)*Fs/K;

maxFreqs=f(index);
delta = QuadraticInterpolation(power_vector,index,f);
ENF_estimation1=(maxFreqs+delta);

% ENF_estimation1=ENF_estimation1/2;

y8 = downsample(ENF,10); % downsample FDR

difLR = numel(y8) - numel(ENF_estimation1);
max_corcoeff=[0 0;0 0];
max_ind=0;
min_std=1;
min_ind=0;
min_Dist=100;
min_ind_Dist=0;
for l=1:difLR
  ENF_temp=0;
  ENF_temp = y8(l:l+1799);
    
R = corrcoef(ENF_temp,ENF_estimation1');    % CORRELATION COEFFICIENT
STD = std(ENF_temp - ENF_estimation1');     %%%% MEAN SQUARED ERROR

if STD < min_std
min_std=STD;
min_ind=l;
end

if (R(1,2))> max_corcoeff(1,2)
max_corcoeff=R;
max_ind=l;
end

end


