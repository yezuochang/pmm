clear all
close all
[f,Y]=ldstone('../demo/HHM1506.s3p');
[f,Y]=ldstone('../demo/ind.s2p');
h=reshape(Y(1,1,:),size(Y,3),1);
semilogy(f,abs(h),'LineWidth',2)
xlabel('freq')
ylabel('Y11')

f=f/max(f)/2/pi;
m=length(f);

omega=2*pi*f(:);
s=j*2*pi*omega;
q=10;
a=[1;zeros(q,1)];
M=[];
for k=1:q
    sk=s.^k;
    M=[M, sk.*h];
end

for k=0:q
    sk=s.^k;
    M=[M, sk];
end

for iter=1:1
    weight = polyval(flipud(a),s);
    M1=M;
    for c=1:size(M1,2)
        M1(:,c)=M1(:,c)./weight;
    end
    Mr=real(M1);
    Mi=imag(M1);
    Hr=real(h./weight);
    Hi=imag(h./weight);
    
    MM=[Mr;Mi];
    HH=[Hr;Hi];
    
    x=MM\HH;
    
    a=x(1:q);
    b=x(q+1:end);
    
    a=[1;a]
    xf=polyval(flipud(b),s)./polyval(flipud(a),s);
    P=M(:,q+1:end)*a;
    Q=M(:,q+1:end)*b;
    xf1=Q./P;

    
end
figure
plot(f,abs(xf), f, abs(h))

P=M(:,q+1:end)*a;
Q=M(:,q+1:end)*b;
M0=M(:,q+1:end);
Jb=M0;
for c=1:size(Jb,2)
    Jb(:,c)=Jb(:,c)./P;
end

b=Jb\h;
figure
plot(f,abs(Jb*b),f,abs(h));


for c=1:100
    P=M0*a;
    Q=M0*b;
    L=-Q./P./P;    
    Ja=M0(:,2:end);
    for d=1:size(Ja,2)
        Ja(:,d)=Ja(:,d).*L;
    end
    J=[Ja Jb];
    x=[a;b];
    xf=Q./P;
    R=xf-h;
    
    A=J.'*J;
    A=A+1*eye(size(A));
    dx=(J.'*R);
    x=x-[0;dx];
end
figure
plot(f,abs(xf),f,abs(h));









