function smax=ss_sigma_max(G,r,p)
if nargin < 3
    p = 1;
end
[n,m]=size(G.B);
s=j*r;
I=eye(n);
H=G.C*((s*I-G.A)\G.B)+G.D;
smax=sort(svd(H),'descend');
smax = smax(p);