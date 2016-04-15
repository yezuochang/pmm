function W=SparseGramian(A,B)
% compute the Gramian W, such that
% AW+WA'=-B*B'
% where A is diagonal
BB=B*B';
[X,Y]=meshgrid(diag(A));

XY=conj(Y)+X;
W=-BB./XY;

