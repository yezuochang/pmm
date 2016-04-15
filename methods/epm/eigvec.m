function [V Unn]=eigvec(A,B,r)
n=length(A);
nv=length(r);
% V=zeros(n,nv);
for c=1:nv
  [L,U]=lu(A-r(c)*B);
  U1=U(1:n-1,1:n-1);
  b1=U(1:n-1,n);
  x1=-U1\b1;
  x=[x1;1];
  V(:,c)=x/norm(x);
  Unn(c,1)=U(n,n);
end
