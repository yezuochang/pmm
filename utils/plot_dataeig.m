function [eigH,detG]=plot_haeig(F,Y,type,color)
if nargin < 3
    type = 'Y';
end
if nargin < 4
    color='red';
end
nf=length(F);
m=size(Y,1);
for c=1:nf
    Y0=Y(:,:,c);
    if strcmp(type,'Y')
        d=eig(Y0+Y0');
    else
        d=eig(Y0'*Y0);
    end
    V(:,c)=d;
end
hold on;
if strcmp(type,'Y')
    semilogx(F,0*F,'k');
else
    semilogx(F,0*F+1,'k');
end

for c=1:m
    semilogx(F,V(c,:),color,'LineWidth',1.5);
end
%axis(limits);

