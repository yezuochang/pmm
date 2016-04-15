function spf_verify(G,W,f,Hd)
nf=length(f);
[n,m]=size(G.B);
Hg=ss_xf(G,f);
Hw=ss_xf(W,f);
if nargin<4
    Hd=[];
end
for c=1:nf
    if G.parametertype == 'Y'
        Hg(:,:,c)=Hg(:,:,c)+Hg(:,:,c)';
    else
        Hg(:,:,c)=eye(m)-Hg(:,:,c)'*Hg(:,:,c);
    end
    Hw(:,:,c)=Hw(:,:,c)'*Hw(:,:,c);    
    if ~isempty(Hd)
        Hd(:,:,c)=Hd(:,:,c)+Hd(:,:,c)';
    end
end
Hg=permute(Hg,[3,1,2]);
Hw=permute(Hw,[3,1,2]);
if ~isempty(Hd)
    Hd=permute(Hd,[3,1,2]);
end
m=min(m,2);
if ~isempty(Hd)
    Hd=abs(real(Hd));
end
Hg=abs(real(Hg));
Hw=abs(real(Hw));
for c=1:m
    for d=1:m
        subplot(m,m,m*(c-1)+d);
%         err=norm(Hg{c,d}./Hg{c,d}-1);
%         fprintf('err = %e\n',err);
        if ~isempty(Hd)
            loglog(f,real(Hd(:,c,d)),'b-',f,real(Hg(:,c,d)),'g-',f,real(Hw(:,c,d)),'r-','LineWidth',2);
            legend('data',sprintf('G+G^H',c,d),sprintf('W^HW',c,d));
        else
            loglog(f,real(Hg(:,c,d)),'g-',f,real(Hw(:,c,d)),'r-','LineWidth',2);
            legend(sprintf('G+G^H',c,d),sprintf('W^HW',c,d));
        end
        maxy=max(max(real(Hg(:,c,d))),max(real(Hw(:,c,d))));
        miny=min(min(real(Hg(:,c,d))),min(real(Hw(:,c,d))));
        if ~isempty(Hd)
            maxy=max(maxy,min(real(Hd(:,c,d))));
            miny=min(miny,min(real(Hd(:,c,d))));
        end
        height=maxy-miny;
        axis([min(f),max(f),miny*0.9,maxy+height*0.2]);
        
        xlabel('Freq');
    end
end