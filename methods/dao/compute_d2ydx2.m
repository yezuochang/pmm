function d2ydx2=compute_d2ydx2(n, m, Jcn)

d2ydx2=compute_d2ydx2_diag(n, m, Jcn);



function d2ydx2=compute_d2ydx2_diag(n, m, Jcn)
nvar=n*m+m*m;
In=speye(n);
Im=speye(m);
Invar=speye(nvar);

parfor c=1:nvar
    if (nvar>1000)
        fprintf('Generating d2ydx2 %d/%d.\n', c, nvar);
    end
    dL  = reshape(Invar(1:n*m,c),m, n);
    dQ  = reshape(Invar(m*n+1:m*n+m*m,c),m,m);
  
    dK1=kronvec(dQ',In);
    dK2=kronvec(Im,dL,1);
    dL1=(kronvec(Im,dQ,1)+kronvec(dQ',Im))/2;
    dN=-(kronvec(dL',In)+kronvec(In,dL,1));
    ix=(c-1)*nvar+1:c*nvar;
    KK=[Jcn*dN+dK1 dK2;sparse(m*m,m*n) dL1]';
    %d2ydx2(ix,:)=KK;
    d2ydx2(c,:)=KK(c,:);
end

function d2ydx2=compute_d2ydx2_full(n, m, Jcn)
nvar=n*m+m*m;
In=speye(n);
Im=speye(m);
Invar=speye(nvar);
d2ydx2={};
parfor c=1:nvar
    if (nvar>1000)
        fprintf('Generating d2ydx2 %d/%d.\n', c, nvar);
    end
    dL  = reshape(Invar(1:n*m,c),m, n);
    dQ  = reshape(Invar(m*n+1:m*n+m*m,c),m,m);
  
    dK1=kronvec(dQ',In);
    dK2=kronvec(Im,dL,1);
    dL1=(kronvec(Im,dQ,1)+kronvec(dQ',Im))/2;
    dN=-(kronvec(dL',In)+kronvec(In,dL,1));
    ix=(c-1)*nvar+1:c*nvar;
    KK=[Jcn*dN+dK1 dK2;sparse(m*m,m*n) dL1]';
    %d2ydx2(ix,:)=KK;
    d2ydx2{c}=KK;
end