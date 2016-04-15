function plotXF(mdl,f,Y,S,supertitle)

H=ss_xf(mdl,f);
if nargin>=4 && ~isempty(S)
    Hs=y2s(H,50);
end
if (nargin<5)
    supertitle='';
end

m=size(H,1);
m=min(m,2);

err=xf_error(Y,H,3);
fprintf('err = %e\n',err);
figure;
H=permute(H,[3,1,2]);
if exist('Y')
    Y=permute(Y,[3,1,2]);
end

for c=1:m
    for d=1:m
        h{c,d}=subplot(m,m,m*(c-1)+d);
        
        if (nargin>2)        
        semilogy(f,abs(Y(:,c,d)),f,abs(H(:,c,d)),f,abs(Y(:,c,d)-H(:,c,d)),'LineWidth',2); %,'-',f,imag(Y(:,c,d)),'-',f,imag(H(:,c,d)),'-');
        %set(gca, 'OuterPosition', [0 0 0.5 .95]) ;
%         maxy=max(max(real(Y{c,d})),max(abs(H(:,c,d))));
%         miny=min(min(real(Y{c,d})),min(abs(H(:,c,d))));
%         
%         height=maxy-miny;
%         axis([min(f),max(f),miny*0.9,maxy+height*0.1]);
        %h_legend=legend('|Y(s): data|','|H(s)|: model','|Y(s)-H(s)|', 'Location','NorthEastOutside');
        %set(h_legend,'FontSize',8);
        else
            plot(f,real(H(:,c,d)),'-',f,imag(H(:,c,d)),'-',f,abs(H(:,c,d)),'-');
        end
        xlabel('Freq');  
        ylabel(sprintf('Y%d%d',c,d));
        %set(h{c,d},'position',[0.4*(d-1) 0.4*(c-1) 0.4 0.4]);
%         subplot(m*m,2,m*(c-1)+d+1);
%         semilogx(f,angle(Y{c,d}),'-',f,angle(H{c,d}));
%         legend('angle(data)','angle(model)');
    end
end
suptitle(sprintf('%s (|E|=%g)',supertitle,err));

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






