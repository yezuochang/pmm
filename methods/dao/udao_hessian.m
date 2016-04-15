function x=udao_hessian(x0,A,B,R,QG,del2,M,parametertype,opts)
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
Invar=speye(nvar);

for c=1:nvar
    dL  = reshape(Invar(1:n*m,c),m, n);
    dQ  = reshape(Invar(m*n+1:m*n+m*m,c),m,m);
  
    dK1=kronvec(dQ',In);
    dK2=kronvec(Im,dL,1);
    dL1=(kronvec(Im,dQ,1)+kronvec(dQ',Im))/2;
    dN=-(kronvec(dL',In)+kronvec(In,dL,1));
    ix=(c-1)*nvar+1:c*nvar;
    KK=[Jcn*dN+dK1 dK2;sparse(m*m,m*n) dL1]';
    %d2ydx2(ix,:)=KK;
    d2ydx2{c}=KK;
end
R2=R'*R;
objective=@(x)udao_objective(x,A,B,R,QG,del2,M,parametertype,Jcn,R2,d2ydx2);

maxiter=optget(opts,'udao_maxiter',500);
if optget(opts,'verbose',0)==1
    displayflag='iter';
else
    displayflag='none';
end
eval_cnt=[];
opt1 = optimset('Jacobian','on','Display',displayflag,'MaxIter',maxiter,...
    'TolX',1e-6*norm(x0),'TolFun',1e-6,'MaxFunEvals',1e5,'LargeScale','on',...
    'GradObj','on','Hessian','on','PrecondBandWidth',0,'DerivativeCheck','off',...
    'Algorithm','interior-point');
[x,res,exitflag,output,grad,hess1] = fminunc(objective,x0,opt1);