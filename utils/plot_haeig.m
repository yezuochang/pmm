function limits=plot_haeig(mdl,limits,step,color)

if nargin<2 || isempty(limits)
    if strcmp(mdl.parametertype,'Y') 
        limits=[0 10 -1 1];
    else
        limits=[0 10 0.8 1.2];
    end
end
if nargin<3 || isempty(step)
    step=1e-2;
end
if nargin<4
    color='red';
end
f2=[];
if nargin < 2 || isempty(limits)
    [r2,f2]=passivity_violation(mdl);
    
    if (length(f2)>0)
        f=[];
        maxf=max(f2)*1.2;
        f3=[0;f2;maxf];
        for c=2:length(f3)
            f=[f f3(c-1)+(0:9)/10*(f3(c)-f3(c-1))];
        end
        limits(1)=0;
        limits(2)=maxf;
    else
        f=limits(1):step:limits(2);
    end   
    
    if length(f) == 0
        r=eig(full(mdl.A));
        maxf=max(imag(r))/2/pi;
        f=maxf*(0:999)/1000;
    end
else
    f=limits(1):step:limits(2);
end

%omg=1e9:1e8:4e11;
omg=2*pi*f;
s=i*omg;
nf=length(s);
[A,B,C,D]=ss_data(mdl);
[n,m]=size(B);
I=eye(n);
for c=1:nf
    H=C*((s(c)*I-A)\B)+D;
    if strcmp(mdl.parametertype,'Y')               
        V(:,c)=eig(H+H');
    elseif strcmp(mdl.parametertype,'S')
        V(:,c)=svd(H);
    else
        assert(0);
    end        
end



hold on;
for c=1:m
    if mdl.parametertype == 'Y'
        plot(f,V(c,:),color,'LineWidth',1.5);
    else
        plot(f,V(c,:),color,'LineWidth',1.5);
    end    
end
if strcmp(mdl.parametertype,'Y')
    plot(f,0*f,'k');
    if ~isempty(f2)
        plot(f2,0*f2,'k*');
    end
else
    plot(f,0*f+1,'k');
    if ~isempty(f2)
        plot(f2,0*f2+1,'k*');
    end
end
axis(limits);
xlabel('Freq (Hz)');
ylabel('\lambda(H+H'')');

