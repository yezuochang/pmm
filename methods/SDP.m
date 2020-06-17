function [G,W,report]=sdp(G0,W0,freq,Yv,opts)

if G0.parametertype == 'S'
    G0 = ss_s2y(G0, 50);
    fprintf('Warning: SDP method only works for hybrid (Y- or Z-) systems.\n');
    fprintf('Will return just AS IS.\n');
    G=G0;
    W=G0;
    return;
end

[R,QG,del2]=preprocess(freq,Yv,G0,opts);
[A,B,C,D]=ss_data(G0);
A=sparse(A);
B=sparse(B);
[n,m]=size(B);
nvar=n*m+m*m;
In=speye(n);
Im=speye(m);
Im2=speye(m^2);

Z=kronvec(A',In)+kronvec(In,A);
[L,U]=lu(Z);
M=-inv(U)*inv(L); % M=-inv(Z)
Jcn=kronvec(B',In)*M;
[C,D]=solve_cvx(n,m,Jcn,R,QG,del2,G0,freq,Yv,opts);
%[C,D]=solve_yalmip(n,m,Jcn,R,QG,del2);
G=G0;
G.C=C;
G.D=D;
if optget(opts,'spf',0); 
    W=spf(G);
else
    W=[];
end

function [C,D]=solve_cvx(n,m,Jcn,R,QG,del2,G,freq,Yv,opts)
if optget(opts,'verbose',0) == 0
    cvx_quiet(true);
else
    cvx_quiet(false);
end
if optget(opts,'enforceDC',0)
    Q = R'*R;
    f = -QG'*R;
    Aeq = [kron(-G.B'*inv(G.A)',eye(m)),eye(m^2)];
    beq = vec(Yv(:,:,1));
    cvx_begin sdp;
        variable Z(n+m,n+m) hermitian;
        Z>=0;
        VC=Jcn*vec(Z(1:n,1:n))+vec(Z(n+1:n+m,1:n));
        VD=vec(Z(n+1:n+m,n+1:n+m))/2;
        y=[VC;VD];
        Aeq*[VC;VD]==beq;
        E=R*y-QG;
        res=E'*E+del2;
        minimize(res);
    cvx_end;
else
    cvx_begin sdp;
        variable Z(n+m,n+m) hermitian;    
        Z-eye(n+m,n+m)*1e-8>=0;
        VC=Jcn*vec(Z(1:n,1:n))+vec(Z(n+1:n+m,1:n));
        VD=vec(Z(n+1:n+m,n+1:n+m))/2;
        y=[VC;VD];
        E=R*y-QG;
        res=E'*E+del2;
        minimize(res);
    cvx_end;
end
C=reshape(VC,m,n);
D=reshape(VD,m,m);

function [C,D]=solve_cvx2(n,m,Jcn,R,QG,del2,G,freq,Yv,opts)
if optget(opts,'verbose',0) == 0
    cvx_quiet(true);
else
    cvx_quiet(false);
end
if optget(opts,'enforceDC',0)
    Q = R'*R;
    f = -QG'*R;
    np = n/m;
    Aeq = [kron(-G.B'*inv(G.A)',eye(m)),eye(m^2)];
    beq = vec(Yv(:,:,1));
    cvx_begin sdp;
        variable Z(n+m,n+m) hermitian;
        Z-eye(n+m,n+m)*1e-6>=0;
        VC=Jcn*vec(Z(1:n,1:n))+vec(Z(n+1:n+m,1:n));
        VD=vec(Z(n+1:n+m,n+1:n+m))/2;
        y=[VC;VD];
        Aeq*[VC;VD]==beq;
        E=R*y-QG;
        res=E'*E+del2;
        minimize(res);
    cvx_end;
else
    np = n/m;
    Z0=zeros(np,np);
    Z1=sparse(n,n);
    for c=1:m
        ix=(c-1)*np+1:c*np;
        Z1(ix,ix)=1;
    end
    [I,J,V]=find(Z1);
    
    cvx_begin sdp;
        variable V(np*np*m,1);
        variable Zb(n,m);
        variable Zd(m,m) hermitian;
        
        Z1=sparse(I,J,V);
        Z=[Z1 Zb;Zb' Zd];
        Z-eye(n+m,n+m)*1e-5>=0;
        VC=Jcn*vec(Z(1:n,1:n))+vec(Z(n+1:n+m,1:n));
        VD=vec(Zd/2);
        y=[VC;VD];
        E=R*y-QG;
        res=E'*E+del2;
        minimize(res);
    cvx_end;
end
C=reshape(VC,m,n);
D=reshape(VD,m,m);

function [C,D]=solve_yalmip(n,m,Jcn,R,QG,del2)
Z=sdpvar(n+m,n+m);
VC=Jcn*vec(Z(1:n,1:n))+vec(Z(n+1:n+m,1:n));
VD=vec(Z(n+1:n+m,n+1:n+m))/2;
y=[VC;VD];
E=R*y-QG;
res=E'*E+del2;
F = [Z >= 0];
solvesdp([],res);

C=reshape(double(VC),m,n);
D=reshape(double(VD),m,m);


