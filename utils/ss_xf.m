function H=ss_xf(mdl,f)
nf=length(f);
A=mdl.A;
B=mdl.B;
C=mdl.C;
D=mdl.D;
if exist('mdl.K')
    K=mdl.K;
else
    K=D*0;
end
%K=mdl.K;
[n,m]=size(B);
I=eye(n);
for c=1:nf
    s=2*pi*i*f(c);
    Hc=D+C*((s*I-A)\B)+s*K;
    H(:,:,c)=Hc;
end
