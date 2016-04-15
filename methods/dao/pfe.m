function G=pfe(W,M)

% Partial Fraction Expansion
% Find G(s) to satisfy
% G(s)+G^H(s)=W^H(s)W(s)
%
A=W.A;
B=W.B;
L=W.C;
Q=W.D;
if nargin < 2
    if max(abs(A-diag(diag(A))))==0
        K=SparseGramian(A, L');
    else
        K=lyap(A',L'*L);
    end
else
    n=size(W.A,1);
    K=reshape(-M*vec(L'*L),n,n);
end
G.A=A;
G.B=W.B;
G.C=B'*K+W.D'*W.C;
G.D=W.D'*W.D/2;
G.parametertype=W.parametertype;

function W=SparseGramian(A,B)
% compute the Gramian W, such that
% AW+WA'=-B*B'
% where A is diagonal
BB=B*B';
[X,Y]=meshgrid(diag(A));

XY=conj(Y)+X;
W=-BB./XY;

