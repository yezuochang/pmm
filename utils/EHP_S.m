function [M,N,J]=EHP_S(S)
A=S.A;
B=S.B;
C=S.C;
D=S.D;
[n,m]=size(B);
I=eye(m);
P=[A 0*A;
   0*A -A'];
L=[B'*0 B'
   C C*0;];
U=[B C'*0;
   B*0 -C'];
Q=[-I D';
   D -I];
M=[P U;L Q];
%M=[A 0*A B; 0*A -A' -C'; C B' D+D'];
N=zeros(2*n+2*m,2*n+2*m);
N(1:n,1:n)=eye(n);
N(n+1:2*n,n+1:2*n)=eye(n);
J=zeros(2*n+m);
J(1:n,n+1:2*n)=-eye(n);
J(n+1:2*n,1:n)=eye(n);
J(2*n+1:2*n+2*m,2*n+1:2*n+2*m)=eye(2*m);

%U=[A 0*A;0*A -A'];
%S=R-[C B']*inv(U)*[B;-C'];
%T=eye(2*n)+[B;-C']*inv(S)*[C B']*inv(U);
%N2=inv(U)*T;