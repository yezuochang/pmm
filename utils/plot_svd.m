function [eigH,detG]=plot_svd(mdl,limits,step,color)

if (nargin<2)
    limits=[1e1 1e12 -0.1 0.1];
end
A=mdl.A;
B=mdl.B;
C=mdl.C;
D=mdl.D;

[n,m]=size(B);


power_f=log10(limits(1)):step:log10(limits(2));
omg=10.^power_f;
%omg=1e9:1e8:4e11;
f=omg;

s=i*omg;
nf=length(s);

I=eye(n);
for c=1:nf
    H=C*((s(c)*I-A)\B)+D;
    svdH(:,c)=svd(H);    
end
hold on;
semilogx(f,0*f,'k');
semilogx(f,1,'-b');
for c=1:m
    semilogx(f,svdH(c,:),color,'LineWidth',1.5);
end
axis(limits);

