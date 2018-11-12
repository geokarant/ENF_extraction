clear all
close all
clc


data_1         = struct2array(load('rec1208-145650_r@441Hz.mat'));        %Load data


ENF         = struct2array(load('ENFdata_FDR.mat'));
 

% signal filtering 

Fs        = 441;                % sampling frequency (Hz)
Ts        = 1/Fs;               % sampling period (seconds)
F_Nyq     = Fs/2;               % Nyquist frequency (Hz)

F_pass1   = 59.95;               %First cut-off frequency (Hz)
F_pass2   = 60.05;               %Second cut-off frequency (Hz)
   
F_pass1_n     = F_pass1/F_Nyq;                                          %Normalized frequency
F_pass2_n     = F_pass2/F_Nyq;                                          %Normalized frequency

N             = 1501;                                                   %Filter order (must be odd number)
win           = hamming(N);                                             %Window function

b             = fir1(N-1,[F_pass1_n F_pass2_n],'bandpass',win,'scale'); %Filter coefficients (impulse response)

data_filtered = filtfilt(b,1,data_1);                                     %Apply filtering


clear b
clear F_Nyq
clear F_pass1
clear F_pass1_n
clear F_pass2
clear F_pass2_n
clear N
clear win

% original signal plot

figure (1)
set(gcf, 'Position', get(0,'Screensize'),'PaperPositionMode','auto')
plot(0:Ts:Ts*(size(data_1,1)-1),data_1)
xlabel('Time (s)','Fontsize',25)
title('Original Signal','Fontsize',25)
set(gca,'Fontsize',25,'Linewidth',2)
grid on

% filtered signal plot

figure (2)
set(gcf, 'Position', get(0,'Screensize'),'PaperPositionMode','auto')
plot(0:Ts:Ts*(size(data_filtered,1)-1),data_filtered,'-r')
xlabel('Time (s)','Fontsize',25)
title('Filtered Signal','Fontsize',25)
set(gca,'Fontsize',25,'Linewidth',2)
grid on


% PSD estimate of both signals

figure (3)
set(gcf, 'Position', get(0,'Screensize'),'PaperPositionMode','auto')
[pxx,f] = pwelch(data_1,500,200,500,Fs);                                % check PSD to see if 60 Hz wins
plot(f,(log10(pxx)))
hold on
[pxx1,f1] = pwelch(data_filtered,500,200,500,Fs);                       % check PSD to see if 60 Hz wins
plot(f1,(log10(pxx1)))
xlabel('Frequency (Hz)','Fontsize',25')
ylabel('Magnitude (dB)','Fontsize',25')
set(gca,'Fontsize',25,'Linewidth',2)
legend('Original Signal','Filtered Signal')
grid on


%%%%%%%%%%%%%%%%   Short Time Fourier Transform  %%%%%%%%%%%%%%%%

%% variable setup
%Fs=441;

frame_len_secs_data1 = 20;                     % frame length equals to 20 seconds
                 
signal_len     = length(data_1);           %signal length 
T              = signal_len/Fs;                 %signal duration


frame_length   = frame_len_secs_data1 * Fs; % frame length in samples


shift_amount   = Fs;                        % 1sec * Fs

  nfft= 4*frame_length;                        % numbers of fft
%nfft = 4*50;
nb_frames = ceil((signal_len - frame_length)/shift_amount)+1;
  
w1 = rectwin(frame_length);                    % Rectangular window 
w2 = rectwin(frame_length/4); 
data_filtered1=data_filtered(1:signal_len);    % dummy variable

%data_filtered1=data_1(1:signal_len);    % dummy variable


%%%%%%%%%%%%%% Manually STFT  %%%%%%%%%%%%%%%%%%%%


% form the stft matrix
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((signal_len-frame_length)/shift_amount);        % calculate the total number of columns
stft = zeros(rown, coln);           % form the stft matrix

% initialize the indexes
indx = 0;
col = 1;
tic
% perform STFT
while indx + frame_length <= signal_len
    
    % windowing
    xw = data_filtered1(indx+1:indx+frame_length).*w1;
    
    
   
    % STFT
      X = fft(xw, nfft);
  
   % Music
   % X =music(xw,2,4);
   
   % rfb
   %X=rfb(xw,1,nfft/4);
   
   % E-sprit
   
  % X = esprit(xw,2,4);
   
   % Capon
  %  X=capon(xw,50,nfft);
  
     % fast Capon

   % X = fast_capon1(y,m,L);
   
   % Daniell method
   % X=daniellse(xw,2,nfft);
   
   % Welch
  
   % X=welchse(xw,w2,1000,nfft);
   
   % Blackman-Tukey method
   
   % X=btse(xw,w2,nfft);
   

    
    % update the stft matrix
      stft(:, col) = X(1:rown);
    
      
    % ESPRIT MATRIX
      %stft(:, col) = X(:,1);
    
    % update the indexes
    indx = indx + shift_amount;
   %tic
    col = col + 1
   %toc 
end

toc
% calculate the time and frequency vectors
t = (frame_length/2:shift_amount:frame_length/2+(coln-1)*shift_amount)/Fs;
f = (0:rown-1)*Fs/nfft;





%% Calculate Power Spectrum %%
for ii=1:col-1;
Power_Spectrum(:,ii)=(abs(stft(:,ii)).^2)/frame_length;
end

power_vector=log10((Power_Spectrum));  % compute power vector

%%%% end manual stft %%%%%%%%%



% %% Mean value for ESPRIT
% 
% FREQS_ESPRIT = stft(:,:).*(Fs/(2*pi));
% FREQS_ESPRIT=abs(FREQS_ESPRIT);
% 
% for jj=1:1800
% ENF_estimation1(1,jj)=(FREQS_ESPRIT(2,jj));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, index] = max(power_vector);       % find indices of max power


delta = QuadraticInterpolation(power_vector,index,f);


maxFreqs=f(index);

ENF_estimation=(maxFreqs+delta);

%ENF_estimation = stft(1,:);
ENF_estimation1=ENF_estimation(1:1800);


tdump=t(1:1800)-9;
% ENF_estimation1=ENF_estimation1/2;      % for scaling 120 hz to 60 hz
% ENF_estimation1=ENF_estimation1/3;      % for scaling 180 hz to 60 hz


%% gia DTW
y8 = downsample(ENF,10);  
dlmwrite('enf.txt',ENF_estimation1,'delimiter','\t','precision',17);
ENF_temp = y8(403:403+1799);
R = corrcoef(ENF_temp,ENF_estimation1')




%%%% cor-coef me ENF ground truth

y8 = downsample(ENF,10);                              % downsample FDR

difLR = numel(y8) - numel(ENF_estimation1);
max_corcoeff=[0 0;0 0];
max_ind=0;
min_std=1;
min_ind=0;

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

%ENF_estimation1=ENF_estimation1+0.05;       %0.05 offset for visualization




ENF_PLOT(:,1)=y8(min_ind:1799+min_ind);



tdump=1:1800;


STFT=STFT-0.05;
WELCH=WELCH+0.05;
Capon_1_180_L1=Capon_1_180_L1+0.1;
ESPRIT=ESPRIT+0.15;

figure (4)
set(gcf, 'Position', get(0,'Screensize'),'PaperPositionMode','auto')
plot(tdump,STFT);
hold on
plot(tdump,ENF_PLOT');
plot(tdump,WELCH);
plot(tdump,Capon_1_180_L1);
plot(tdump,ESPRIT);
xlabel('Time (s)','Fontsize',25)
ylabel('Frequency (Hz)','Fontsize',25)
%xlim([0 1800])
ylim([59.8 60.2])
legend('STFT','Ground truth','Welch','Capon ','ESPRIT')
%title('ENF Estimation Data 1 - 120 Hz - ESPRIT','Fontsize',25)
set(gca,'YDir','Reverse','Fontsize',25,'Linewidth',2)
grid on
%export_fig newFINAL1.png
export_fig newFINAL.png -transparent





