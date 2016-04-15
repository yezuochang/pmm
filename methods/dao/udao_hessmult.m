function x=udao_hessmult(x0,A,B,R,QG,del2,M,parametertype,opts)
global eval_cnt
if opts.verbose 
    fprintf('Preparing for optimization.\n');
end
[n,m]=size(B);
nvar=n*m+m*m;
In=speye(n);
Im=speye(m);
Im2=speye(m^2);
Jcn=kronvec(B',In)*M;

d2ydx2=compute_d2ydx2(n, m, Jcn);

R2=R'*R;
objective=@(x)udao_objective_hessmult(x,A,B,R,QG,del2,M,parametertype,Jcn,R2,d2ydx2);

maxiter=optget(opts,'udao_maxiter',500);
if optget(opts,'verbose',0)==1
    displayflag='iter';
else
    displayflag='none';
end
MaxPCGIter=round(min(nvar/2, 1000));
eval_cnt=[];
opt1 = optimset('Jacobian','on','Display',displayflag,'MaxIter',maxiter,...
    'MaxPCGIter',MaxPCGIter,...
    'TolX',1e-6*norm(x0),'TolFun',1e-6,'MaxFunEvals',1e5,'LargeScale','on',...
    'GradObj','on','Hessian','on','HessMult',@hessmult,'Preconditioner',@precon,'PrecondBandWidth',0,'DerivativeCheck','off',...
    'Algorithm','interior-point');
[x,res,exitflag,output,grad,hess1] = fminunc(objective,x0,opt1);

function W=hessmult(Hinfo,Y)
R=Hinfo.R;
dydx=Hinfo.dydx;
W1=dydx*Y;
W1=R*W1;
W1=R'*W1;
W1=dydx'*W1;
W1=2*W1;
%W1=2*T'*(dxdy*Y);
W2=Hinfo.hess2*Y;
W=W1+W2;
%W=(Hinfo.hess1+Hinfo.hess2)*Y;

function [R,pvec] = precon(Hinfo,varargin)
%H=Hinfo.hess1+Hinfo.hess2;
n=size(Hinfo.hess2,1);
epsi=1e-16;
%T=Hinfo.T;
R=Hinfo.R;
dydx=Hinfo.dydx;

dnrms=sum(R.*R).*sum(dydx.*dydx);
%dnrms = sum(T.*T);
%dnrms = sqrt(sum(H.*H))';
d = max(sqrt(dnrms),epsi);
R = sparse(1:n,1:n,full(d));
pvec = (1:n);

