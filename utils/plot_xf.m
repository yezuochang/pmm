function plotXF(G,F,H)

H1=ss_xf(G,F);
if nargin<3
    H=[];
end

m=size(H,1);
m=min(m,2);

if exist('H')
    err1=xf_error(H1,H,2);
%     fprintf('err1 = %e\n',err1);
    H=permute(H,[3,1,2]);
end
H1=permute(H1,[3,1,2]);

for c=1:m
    for d=1:m
        subplot(m,m,m*(c-1)+d);
        
        if (nargin>2)        
            h=abs(H(:,c,d));
            h1=abs(H1(:,c,d));
            e1=abs(H1(:,c,d)-H(:,c,d));
            plot(F,log10(h),...
                 F,log10(h1),...
                 F,log10(e1),...
                 'LineWidth',1.5);
        else
            plot(f,real(H(:,c,d)),'-',f,imag(H(:,c,d)),'-');
        end
        miny=min(e1)*0.5;
        maxy=max(h)*2;
        miny=log10(miny);
        maxy=log10(maxy);
        minF=min(F);
        maxF=max(F);
        xlabel('Freq');  
        if (G.parametertype=='Y')
            ylabel(sprintf('log_{10}|Y_{%d%d}|',c,d));
        else
            ylabel(sprintf('log_{10}|S_{%d%d}|',c,d));
        end
        axis([minF maxF miny maxy]); 
    end
end
%suptitle(sprintf('%s (|E|=%g)',supertitle,err));
% subtitle(' ');

return;

if (nargin>=4 && ~isempty(S))
figure
Hs=permute(Hs,[3,1,2]);
S=permute(S,[3,1,2]);
for c=1:m
    for d=1:m
        subplot(m,m,m*(c-1)+d);
        
        if (nargin>2)
        err=norm(S(:,c,d)./Hs(:,c,d)-1);
        fprintf('err = %e\n',err);
        plot(f,abs(S(:,c,d)),'-',f,abs(Hs(:,c,d)));
%         maxy=max(max(abs(S(:,c,d))),max(abs(Hs(:,c,d))));
%         miny=min(min(abs(S(:,c,d))),min(abs(Hs(:,c,d))));
%         height=maxy-miny;
%         axis([min(f),max(f),miny*0.9,maxy+height*0.1]);
%         legend('abs(data)','abs(model)');
        else
            loglog(f,abs(Hs(:,c,d)))
        end
        xlabel('Freq');        
%         subplot(m*m,2,m*(c-1)+d+1);
%         semilogx(f,angle(Y{c,d}),'-',f,angle(H{c,d}));
%         legend('angle(data)','angle(model)');
    end
end
end






