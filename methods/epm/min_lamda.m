function lamda=min_lamda(omg)
global mmdl;
A=mmdl.A;
B=mmdl.B;
C=mmdl.C;
D=mmdl.D;

[n,m]=size(B);

s=i*omg;

I=eye(n);
H=C*((s*I-A)\B)+D;
G=H+H';
lamda=min(eig(G));