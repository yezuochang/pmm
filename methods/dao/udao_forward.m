function W=udao_forward(G)
if G.parametertype=='Y'
    W=spf(G);
else
    H=pfe(G);
    H.C=-H.C;
    H.D=1/2*eye(size(H.D))-H.D;
    W=spf(H);
end