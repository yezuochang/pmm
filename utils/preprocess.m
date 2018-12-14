function [R,QG,del2]=preprocess(freq,Yv,G0,opts)

if ~isreal(G0.A)
    G0=ss_real(G0);
end
wgt_scheme = optget(opts,'wgt_scheme',3);

[n,m]=size(G0.B);
m2=m*m;
nvar=n*m+m2;
In=speye(n);
Im=speye(m);
Im2=speye(m^2);
nf=length(freq);

[A,B,C,D]=ss_data(G0);
A=sparse(A);
B=sparse(B);
%F=repmat(0,m2*nf,nvar);
%F=zeros(m2*nf,nvar);
T=[];
G=zeros(m2*nf,1);

for k=1:nf
    s=2*pi*j*freq(k);    
    JJ=[kronvec(Im,(s*In-A)\B) kronvec(Im,Im)];
    switch wgt_scheme
        case 1
            wgt=ones(m2,1);
        case 2
            wgt=1./vec(abs(Yv(:,:,k)));            
        case 3
            wgt=1./vec(sqrt(abs(Yv(:,:,k))));
        case 4
            wgt=ones(m2,1)./norm(Yv(:,:,k));
        case 5
            wgt=ones(m2,1)./sqrt(norm(Yv(:,:,k)));
        otherwise
            wgt=ones(m2,1);
    end
    %wgt=min(wgt,100);
    wgt=spdiags(wgt, 0, m2, m2);

    Fk=wgt*JJ;
    Gk=wgt*vec(Yv(:,:,k));
    [I1,J1,V1]=find(Fk);   
    if k == 1
        I=repmat(I1(:),1,nf);
        J=repmat(J1(:),1,nf);
        V=repmat(V1(:),1,nf);        
    else
        I(:,k)=I1+(k-1)*m2;
        J(:,k)=J1;
        V(:,k)=V1;
    end
    %F=[F;Fk];
    %F((k-1)*m2+1:k*m2,1:nvar)=Fk;
    G((k-1)*m2+1:k*m2)=Gk;
end

F=sparse(I,J,V,m2*nf,nvar);
F=[real(F);imag(F)];
G=[real(G);imag(G)];

[QG,R]=qr(F,G,0);
%QG=Q'*G;
del2=G'*G-QG'*QG;
R=sparse(R);