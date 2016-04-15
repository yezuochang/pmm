function [G,W,report]=cdao(G0,W0,freq,Yv,opts)
if ~optget(opts,'enable_udao',1)
    G=G0;
    W=[];
    report.time=0;
    return;
end

if ~ss_isreal(G0)
    G0=ss_complex(G0);
    G0=ss_real(G0);
    G0.C=G0.C*0.9999;
    G0.D=G0.D*0.9999;
end

[R,QG,del2]=preprocess(freq,Yv,G0,opts);  % used to be preprocess2 ??? 
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
M=inv(U)*inv(L); % M=inv(Z)
Jcn=kronvec(B',In)*M;
R2=R'*R;

[x0,info]=cdao_initialvalue(G0,Jcn);
if info < 0
    fprintf('Warning: CDAO cannot get an initial solution.\n');
    fprintf('Will abort without changing the system.\n');
    G=G0;
    W=[];
    return;
end

objective=@(x)cdao_objective(x,R,QG,del2,R2);
hess=[2*R2 R2*0;R2*0 R2*0];
hessfunc=@(x,lambda)hess;
    
if G0.parametertype == 'Y'   
    constr=@(x)cdao_Y_constr(x,A,B,M,Jcn,n,m,G0.parametertype);
else
    constr=@(x)cdao_S_constr(x,A,B,M,Jcn,n,m,G0.parametertype);
end
maxiter=100;
opt1 = optimset('Display','iter','MaxIter',maxiter,...
    'TolX',1e-5,'TolFun',1e-3,'MaxFunEvals',1e5,'LargeScale','off',...
    'GradObj','on','Hessian','on','HessFcn',hessfunc,'GradConstr','on',...
    'TolCon',1e-5,'DerivativeCheck','off', 'Algorithm','interior-point');

[x,fval] = fmincon(objective,x0,[],[],[],[],[],[],constr,opt1);


x1=x(1:nvar);
x2=x(nvar+1:2*nvar);
G=G0;
G.C=reshape(x1(1:n*m),m, n);
G.D=reshape(x1(m*n+1:m*n+m*m),m,m);
W=G0;
W.C=reshape(x2(1:n*m),m, n);
W.D=reshape(x2(m*n+1:m*n+m*m),m,m);

function [x0,info]=cdao_initialvalue(G,Jcn)

if G.parametertype=='Y'
    W=spf(G);
    x0=[vec(G.C);vec(G.D);vec(W.C);vec(W.D)];
    info=1;
else
    [x0,fval]=cdao_forward_S(G,Jcn);
    x0=[vec(G.C);vec(G.D);x0];
    if abs(norm(fval))>1
        info=-1;
    else
        info=1;
    end
end


function [res,grad]=cdao_objective(x,R,QG,del2,R2)
nvar=length(x)/2;
x1=x(1:nvar);
E=R*x1-QG;
res=E'*E+del2;
%grad=2*x1'*R2-2*QG'*R;
grad=2*E'*R;
grad=[grad zeros(size(grad))];

function [c,ceq,gradc,gradceq]=cdao_Y_constr(x,A,B,M,Jcn,n,m,parametertype)
nvar=n*m+m*m;
Invar=speye(nvar);
In=speye(n);
Im=speye(m);
c=0;
gradc=zeros(2*nvar,1);
x1=x(1:nvar);
x2=x(nvar+1:2*nvar);
W.A=A;
W.B=B;
W.C=reshape(x2(1:n*m),m, n);
W.D=reshape(x2(m*n+1:m*n+m*m),m,m);
W.parametertype=parametertype;
G=udao_backward(W,M);
ceq=x1-[vec(G.C);vec(G.D)];
gradceq=[Invar -grad_pfe(W,Jcn)]';
return;
%ceq=delx'*delx;

L=W.C;
Q=W.D;
[m,n]=size(W.C);
In=speye(n);
Im=speye(m);

N=-(kronvec(L',In)+kronvec(In,L,1));

K1=kronvec(Q',In);
K2=kronvec(Im,L,1);
L1=(kronvec(Im,Q,1)+kronvec(Q',Im))/2;

dydx=[Jcn*N+K1 K2;sparse(m*m,m*n) L1];
gradceq=[Invar -dydx]';
%gradceq=2*delx'*[Invar -dydx];
%gradceq=gradceq';
%gradceq=2*delx'*[Invar -dydx];

function [c,ceq,gradc,gradceq]=cdao_S_constr(x,A,B,M,Jcn,n,m,parametertype)
nvar=n*m+m*m;
Invar=speye(nvar);
In=speye(n);
Im=speye(m);
c=0;
gradc=zeros(2*nvar,1);
x1=x(1:nvar);
x2=x(nvar+1:2*nvar);

W1=ssm(A,B,reshape(x1(1:n*m),m, n),reshape(x1(m*n+1:m*n+m*m),m,m),'S');
W2=ssm(A,B,reshape(x2(1:n*m),m, n),reshape(x2(m*n+1:m*n+m*m),m,m),'S');
G1=pfe(W1);
G2=pfe(W2);
y1=[vec(G1.C);vec(G1.D)];
y2=[vec(G2.C);vec(G2.D)];
b=[zeros(n*m,1);vec(eye(m))/2];
ceq=y1+y2-b;
grad1=grad_pfe(W1,Jcn);
grad2=grad_pfe(W2,Jcn);
gradceq=[grad1 grad2]';

function dydx=grad_pfe(W,Jcn)
[n,m]=size(W.B);
L=W.C;
Q=W.D;
[m,n]=size(W.C);
In=speye(n);
Im=speye(m);

N=-(kronvec(L',In)+kronvec(In,L,1));

K1=kronvec(Q',In);
K2=kronvec(Im,L,1);
L1=(kronvec(Im,Q,1)+kronvec(Q',Im))/2;

dydx=[Jcn*N+K1 K2;sparse(m*m,m*n) L1];

function [x,fval]=cdao_forward_S(G,Jcn)
H1=pfe(G);
H2=H1;
H2.C=-H1.C;
H2.D=1/2*eye(size(H1.D))-H1.D;
W=spf(H2);
x0=[vec(W.C);vec(W.D)];

maxiter=20;
opt1 = optimset('Display','iter','MaxIter',maxiter,...
    'TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',1e5,...
    'Jacobian','on','DerivativeCheck','off');
objective=@(x)cdao_forward_S_objective(x,H1,Jcn); 
[x,fval] = fsolve(objective,x0,opt1);


function [fval,jacob]=cdao_forward_S_objective(x,H1,Jcn)
[n,m]=size(H1.B);
In=speye(n);
Im=speye(m);

W=ssm(H1.A,H1.B,reshape(x(1:n*m),m, n),reshape(x(m*n+1:m*n+m*m),m,m),'S');
H2=pfe(W);
y1=[vec(H1.C);vec(H1.D)];
y2=[vec(H2.C);vec(H2.D)];
b=[zeros(n*m,1);vec(eye(m))/2];
fval=y1+y2-b;
jacob=full(grad_pfe(W,Jcn));