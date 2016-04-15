function [W,report]=spf(G)

G=ss_real(G);
assert(ss_isreal(G));

[W,K,report]=spf_sub(G,0);
if report >= 0
    return;
end

[n,m]=size(G.B);
G1=G;
r=1e-5;
alpha=10;
done = 0;
while ~done  
    r=r*alpha;
    G1.D=G.D+eye(m)*r;
    [W,K0,report]=spf_sub(G1,0);
    if report >= 0
        done=1;
    end    
end
[W,report]=spf_sub(G,1,K0);

function [W,K,report]=spf_sub(G,do_newton,K0)


A=G.A;
B=G.B;
C=G.C;
D=G.D;
parametertype=G.parametertype;


R=D+D';
Rinv=inv(R);
Q=chol(R);
Qinv=inv(Q);
Ah=B*Rinv*C-A;
Bh=B*Qinv;
Ch=-C'*Rinv*C;
Ch=(Ch+Ch')/2;
if do_newton
    [K,report]=care_newton(Ah,Bh,Ch,[],[],[],K0); % iterative method as an alternative
    L=Qinv'*C-Qinv'*B'*K;
    W=ssm(A,B,L,Q,parametertype);
else
    [K,L1,G1,report]=care(Ah,Bh,Ch);
    if report >= 0
        L=Qinv'*C-Qinv'*B'*K;
        W=ssmodel(A,B,L,Q,parametertype);
    else
        W=[];
    end
end    





