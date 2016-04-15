function G=udao_backward(W,M)
if W.parametertype=='Y'
    G=pfe(W,M);
else
    H=pfe(W);
    H.C=-H.C;
    H.D=1/2*eye(size(H.D))-H.D;
    G=spf(H);
end