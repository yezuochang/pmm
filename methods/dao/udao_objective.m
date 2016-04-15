function [res,J,hess]=udao_objective(x,A,B,R,QG,del2,M,parametertype,Jcn,R2,d2ydx2)
global  hess_old eval_cnt
[n,m]=size(B);
nm=n*m;
m2=m*m;
nvar=n*m+m*m;
W.A=A;
W.B=B;
W.C=reshape(x(1:n*m),m, n);
W.D=reshape(x(m*n+1:m*n+m*m),m,m);
W.parametertype=parametertype;
G=pfe(W,M);

x1=[vec(G.C);vec(G.D)];
E=R*x1-QG;
res=E'*E+del2;
if nargout == 1
    return;
end
L=W.C;
Q=W.D;
[m,n]=size(W.C);
In=speye(n);
Im=speye(m);

%N=-(kronvec(L',In)+kronvec(In,L,1));

K1=kronvec(Q',In);
K2=kronvec(Im,L,1);
L1=(kronvec(Im,Q,1)+kronvec(Q',Im))/2;

dfdy=2*E'*R;
% dydx=[Jcn*N+K1 K2;sparse(m*m,m*n) L1];
% J=dfdy*dydx;
v1=dfdy(1:nm);
v2=dfdy(nm+1:end);


%J=[(v1*Jcn)*N+v1*K1 v1*K2+v2*L1];


Lv1=reshape(v1*Jcn,n,n);
y1=-vec(L*Lv1);
y2=-vec(Lv1*L');
P=reshape(1:n*m,n,m);
P=reshape(P',n*m,1);
y3=y1+y2(P);
J=[y3'+v1*K1 v1*K2+v2*L1];

if nargout == 2
    return;
end
% if ~isempty(eval_cnt) && eval_cnt < 10
%     eval_cnt=eval_cnt+1;
%     hess=hess_old;
%     return;
% end
N=-(kronvec(L',In)+kronvec(In,L,1));
dydx=[Jcn*N+K1 K2;sparse(m*m,m*n) L1];
dydx=full(dydx);
T=R*dydx;
hess1=2*T'*T;
hess2=sparse(nvar,nvar);
for c=1:nvar
    hess2(:,c)=d2ydx2{c}*dfdy';
end
hess=hess1+hess2;

% maxhess=max(max(abs(hess)));
% ix=find(abs(hess+eye(size(hess))*maxhess)<maxhess*1e-5);
% hess(ix)=0;
% hess=sparse(hess);
% fprintf('nnz(hess)/n^2 = %e\n', nnz(hess)/prod(size(hess)));
hess_old=hess;
eval_cnt=0;
