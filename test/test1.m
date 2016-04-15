clear all
close all
[f,Y]=ldstone('../demo/HHM1506.s3p');
% [f,Y]=ldstone('../demo/ind.s2p');
h=reshape(Y(1,1,:),size(Y,3),1);
semilogy(f,abs(h),'LineWidth',2)
xlabel('freq')
ylabel('Y11')

f=f/max(f)/2/pi;
m=length(f);

omega=2*pi*f(:);
s=j*2*pi*omega;
q=12;

p=-[0.01:0.01:1]';
p=[p/100+j*p;p/100-j*p];
q=length(p);

for c=1:10
    M=[];
    M1=[];
    for c=1:length(p)
        M(:,c)=1./(s-p(c));
        M1(:,c)=-M(:,c).*h;
    end
    A=[M,M1];
    x=A\h;
    r=x(q+1:end);

    p=eig(diag(p)-ones(q,1)*r');
    p
end

M=[];
for c=1:length(p)
    M(:,c)=1./(s-p(c));
end
r=M\h;
xf=M*r;

plot(f,abs(xf),f,abs(h))









