function [R,QG,del2]=preprocess(freq,Yv,G0)

if ~isreal(G0.A)
    G0=ss_real(G0);
end
wgt_scheme = 3;

[n,m]=size(G0.B);
m2=m*m;
In=speye(n);
Im=speye(m);
Im2=speye(m^2);
nf=length(freq);

[A,B,C,D]=ss_data(G0);
A=sparse(A);
tic;
for k=1:nf
    s=2*pi*j*freq(k);    
    J=[kronvec(Im,(s*In-A)\B) kronvec(Im,Im)];
    %J=full(J);
    if (wgt_scheme == 0)
        wr=ones(m2,1);
        wi=ones(m2,1);
    elseif wgt_scheme == 1
        wr=1./vec(sqrt(abs(real(Yv(:,:,k)))));
        wi=1./vec(sqrt(abs(imag(Yv(:,:,k)))));
    elseif wgt_scheme == 2
        wr=1./vec(abs(real(Yv(:,:,k))));
        wi=1./vec(abs(imag(Yv(:,:,k))));    
    elseif wgt_scheme == 3
        wr=1./vec(sqrt(abs(Yv(:,:,k))));
        wi=1./vec(sqrt(abs(Yv(:,:,k))));        
    elseif wgt_scheme == 4
        wr=1./vec(abs(Yv(:,:,k)));
        wi=1./vec(abs(Yv(:,:,k)));
    end
    wr=min(wr,100);
    wi=min(wi,100);
    wr=spdiags(wr, 0, m2, m2);
    wi=spdiags(wi, 0, m2, m2);
    %wi=wi/10000;
    F(:,:,k)=full(wr*real(J));
    F(:,:,k+nf)=full(wi*imag(J));    
    G(:,k)=wr*vec(real(Yv(:,:,k)));
    G(:,k+nf)=wi*vec(imag(Yv(:,:,k)));
end
toc
F=permute(F,[3 2 1]);
G=permute(G,[2,1]);
nvar=n*m+m*m;
del2=0;
FF=[];
GG=[];
for c=1:m2
    FF=[FF;F(:,:,c)];
    GG=[GG;G(:,c)];
end
%    ix=(c-1)*nvar+1:c*nvar;
[Q,R]=qr(FF,0);
QG=Q'*GG;
del2=GG'*GG-QG'*QG;
%end
R=sparse(R);