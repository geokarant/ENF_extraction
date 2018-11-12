function phi=fast_capon1(y,m,L)
%
% The Capon spectral estimator.
%
% phi=capon(y,m,L);
%
%    y   <- the data vector (length N)
%    m   <- the length of the Capon filter
%    L   <- the number of estimated spectral samples
%    phi -> the estimated spectrum
%
% Copyright 2018 by C. Kotropoulos

y=y(:);
N=length(y);       % data length
M=m+1;
sigma2=1e-8;
%% form the first column of the sample covariance matrix
for ii=1:M
     r_M(ii,1)=0;
     for t = M:1:N,
        r_M(ii,1)=r_M(ii,1)+y(t-ii+1,1)*conj(y(t,1));
     end
end
r_M=r_M/(N-m);
r_M(1,1)=r_M(1,1)+sigma2;
%% compute the inverse of R using GS
[m1,m2,m3,m4]=topinv(r_M,r_M');
%Generators of Rinv
t_M=m1{1,1}; %column vector
w_Mminus1=t_M(2:end);
mm2=m2{1,2}; 
t_M=t_M*sqrt(mm2(1,1));
s_M=vertcat(0,conj(flipud(w_Mminus1)))*sqrt(mm2(1,1));
I_Mminus1=eye(M-1);
Z_M=zeros(M,M);
Z_M(2:end,1:end-1) = I_Mminus1; % unit lower shift matrix 
Z_M2 = Z_M;
Z_M2(1,M) = 1;
I=eye(M);
J = fliplr(I); % exchange matrix
% gg1=gallery('krylov',Z_M,t_M,M);
% gg2=gallery('krylov',Z_M,s_M,M);
% IR=gg1*gg1'-gg2*gg2';
% % compute the spectrum
% phi=zeros(L,1);
% for k = 1 : L, 
%    a=exp(- j*2*pi*(k-1)/L*[0:m].');      % form the a(w) vector
%    phi(k)=real(a'*IR*a);
% end
% phi=M./phi;
%% Computation of the polynomial in the denominator of IAA
vv=1:1:M;
tt=flipud(t_M);
tt_M=vv'.* tt; %\tilde{t_M}
ss=flipud(s_M);
ss_M=vv'.* ss; %\tilde{s_M}
cc1=gallery('krylov',Z_M,tt_M,M)*conj(t_M);
cc2=gallery('krylov',Z_M,ss_M,M)*conj(s_M);
cc_=conj((cc1-cc2)); %\underline{c}; 
%cc_=[cc_{-M+1}, \ldots, cc_{-1}, cc_{0}]^T
%cc=vertcat(cc_(M), flipud(cc_(1:M-1)), zeros(L-2*M+1,1), flipud(conj(cc_(1:M-1))));
cc=vertcat(cc_(M), flipud(conj(cc_(1:M-1))),zeros(L-2*M+1,1),cc_(1:M-1));
denominator=fft(cc,L);
% compute the spectrum
phi=M./denominator;
