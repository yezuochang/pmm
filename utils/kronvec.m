function M=kronvec(A,B,transpose)
% vec(AXB)= M vec(X)
% transpose=0: vec(A*X*B) = kron(B',A)*vec(X)
% transpose=1: vec(A*X'*B) = P kron(A,B')*vec(X)
% see http://en.wikipedia.org/wiki/Kronecker_product
if nargin < 3
    transpose=0;
end
if transpose==0
    M=kron(B.',A);
else
    M=kron(A.',B).';
    n=size(A,1);
    m=size(B,2);
    P=reshape(1:n*m,m,n);
    P=reshape(P',n*m,1);
    M=M';
    M=M(:,P);
    M=M';
    %M=M(P,:);
end
%M=sparse(M);