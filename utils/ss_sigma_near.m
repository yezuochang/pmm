function [smax ,num]=ss_sigma_near(G,r,tol)
if nargin < 3
    tol = 0;
end
[n,m]=size(G.B);
s=j*r;
I=eye(n);
H=G.C*((s*I-G.A)\G.B)+G.D;
smax=svd(H) - 1 + tol;

if 0 == tol
    abstol = 1e-6;
else
    abstol = min(0.1*tol,1e-6);
end
ix = find(abs(smax) < abstol);
num = max([length(ix),1]);

[C,IX] = min(abs(smax));
smax = smax(IX);
% smax = abs(smax(IX));