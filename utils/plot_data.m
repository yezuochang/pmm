function plotXF(f,Y)

m=size(Y,1);
m=min(m,2);

Y=permute(Y,[3,1,2]);


for c=1:m
    for d=1:m
        h{c,d}=subplot(m,m,m*(c-1)+d);
        semilogy(f,abs(Y(:,c,d)),'LineWidth',2); %,'-',f,imag(Y(:,c,d)),'-',f,imag(H(:,c,d)),'-');        
    end
end
