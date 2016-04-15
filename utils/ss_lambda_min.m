function lmin=ss_lambda_min(G,r,p)
if nargin < 3
    p = 1;
end
[n,m]=size(G.B);
s=j*r;
I=eye(n);
H=G.C*((s*I-G.A)\G.B)+G.D;
lmin=sort(real(eig(H+H')));
lmin = lmin(p);