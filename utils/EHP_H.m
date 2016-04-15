function [M,N,J]=EHP_H(mod)
A=mod.A;
B=mod.B;
C=mod.C;
D=mod.D;
[n,m]=size(B);
I=eye(m);
M=[A 0*A B; 0*A -A.' -C.'; C B.' D+D.'];
%M=[A 0*A B; 0*A -A' -C'; C B' D+D'];
N=zeros(2*n+m,2*n+m);
N(1:n,1:n)=eye(n);
N(n+1:2*n,n+1:2*n)=eye(n);
J=zeros(2*n+m);
J(1:n,n+1:2*n)=-eye(n);
J(n+1:2*n,1:n)=eye(n);
J(2*n+1:2*n+m,2*n+1:2*n+m)=eye(m);

%U=[A 0*A;0*A -A'];
%S=R-[C B']*inv(U)*[B;-C'];
%T=eye(2*n)+[B;-C']*inv(S)*[C B']*inv(U);
%N2=inv(U)*T;

