function [X,report]=care_newton(A,B,Q,R,S,E,X)

% Petter Benner, etc, An Exact Line Search Method for Solving Generalized 
% Continuous-Time Algebraic Riccati Equations, 1998, IEEE TAC

[n,m]=size(B);

if nargin < 4 || isempty(R)
    R=eye(m);
end
if nargin < 5 || isempty(S)
    S=zeros(n,m);    
end
if nargin < 6 || isempty(E)
    E=eye(n);
end
if nargin < 7 || isempty(X)
    %X=care(A,B,Q,R,S,E);
    X=eye(n);
end
Rinv=inv(R);
abs_err=1;
rel_err=1;
tol=1e-13;
rel_tol=1e-5;
maxiter=1000;
k=0;
while abs_err>tol && rel_err > rel_tol && k < maxiter
    K=Rinv*(B'*X*E+S');
    Z=care_eval(A,B,Q,R,S,E,X);
    N=lyap((A-B*K)',Z);
    X=line_search(A,B,Q,R,S,E,X,N);
    err0=abs_err;
    abs_err=norm(Z);
    rel_err=abs(1-abs_err/err0);
    fprintf('abs_err = %e, rel_err = %e\n', abs_err, rel_err);
    k=k+1;
end
if abs_err > tol
    fprintf('Warning: Newton''s method for solving CARE equation stops when \n');
    fprintf('abusolute error = %e. Result could be inaccurate.\n', abs_err);
    report=-1;
else
    report=1;
end

function X=line_search(A,B,Q,R,S,E,X,N)
min_err=1e10;
best_lambda=1;
for c=1:100
    lambda=(c-1)/100;
    Z=care_eval(A,B,Q,R,S,E,X+lambda*N);
    r(c)=norm(Z);
    if (r(c)<min_err)
        min_err=r(c);
        best_lambda=lambda;
    end
end
% c=1:100;
% plot(c,r);
% pause
X=X+best_lambda*N;
%fprintf('lambda = %e, min_err = %e\n', best_lambda,min_err);

function Z=care_eval(A,B,Q,R,S,E,X)
[n,m]=size(B);
if nargin <5 || isempty(R)
    R=eye(m);
end
if nargin < 6  || isempty(S)   
    S=zeros(n,m);    
end
if nargin < 7 || isempty(E)    
    E=eye(n);
end

T=E'*X*B+S;
Rinv=inv(R);
Z=A'*X*E+E'*X*A-T*Rinv*T'+Q;
Z=(Z+Z')/2;