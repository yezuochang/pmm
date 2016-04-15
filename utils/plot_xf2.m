function plotXF2(G1,G2,F,H)
F=F(1:3:end);
H=H(:,:,1:3:end);
H1=ss_xf(G1,F);
H2=ss_xf(G2,F);
if nargin<3
    H=[];
end

m=size(H,1);
m=min(m,2);




if exist('H')
    err1=xf_error(H1,H,3);
    err2=xf_error(H2,H,3);
    fprintf('err1 = %e, err2 = %e\n',err1, err2);
    H=permute(H,[3,1,2]);
end
H1=permute(H1,[3,1,2]);
H2=permute(H2,[3,1,2]);

figure;
for c=1:m
    for d=1:m
        %h{c,d}=
        subplot(m,m,m*(c-1)+d);
        
        if (nargin>2)
            h=abs(H(:,c,d));
            h1=abs(H1(:,c,d));
            h2=abs(H2(:,c,d));
            e1=abs(H1(:,c,d)-H(:,c,d));
            e2=abs(H2(:,c,d)-H(:,c,d));
            plot(F,log10(h),...
                 F,log10(h1),...
                 F,log10(h2),...
                 F,log10(e1),...
                 F,log10(e2),':',...
                 'LineWidth',1.5); 
        else
            plot(f,real(H(:,c,d)),'-',f,imag(H(:,c,d)),'-');
        end
        miny=min(min(e1),min(e2))*0.5;
        maxy=max(h)*2;
        miny=log10(miny);
        maxy=log10(maxy);
        minF=min(F);
        maxF=max(F);
        xlabel('Freq');  
        if (G1.parametertype=='Y')
            ylabel(sprintf('log_{10}|Y_{%d%d}|',c,d));
        else
            ylabel(sprintf('log_{10}|S_{%d%d}|',c,d));
        end
        axis([minF maxF miny maxy]);        
    end
end


