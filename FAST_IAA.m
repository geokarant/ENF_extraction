%% FIAA
value=20e10;
threshold=1e-6;
sigma2=1e-8;
M=8;    
K=2*M;
% initialization
xw = x(1:M);
iaa_initial_est = ifft(xw,K)*(K/M);
p=abs(iaa_initial_est).^2; %K X 1
iter=1;
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
       disp(['iter:  ', num2str(iter), ' criterion: ', num2str(value)]);
       p=newp;
       iter=iter+1;
end   