function [G,W]=udao(G0,W0,freq,Yv,opts);
if G0.parametertype ~= 'Y';
    fprintf('Warning: UDAO works only for hybrid systems.\n');
    G=ss_real(G0);
    W=[];
    report.time=0;
    return;
end

if ~ss_isreal(G0)
    G0=ss_real(G0);
end

r2=passivity_violation(G0);
if ~isempty(r2)
    figure;
    plot_haeig(G0);
    error('Error: need to fix the passivity first.');    
end

[R,QG,del2]=preprocess(freq,Yv,G0,opts);

[A,B,C,D]=ss_data(G0);
A=sparse(A);
B=sparse(B);
[n,m]=size(B);
nvar=n*m+m*m;
In=speye(n);
Im=speye(m);
Im2=speye(m^2);

Z=kronvec(A',In)+kronvec(In,A);
[L,U]=lu(Z);
M=inv(U)*inv(L); % M=inv(Z)

W=spf(G0);
x0=[vec(W.C);vec(W.D)];

opt_method=optget(opts,'udao_opt_method','hess');

if strcmp(opt_method,'hess')
    x=udao_hessian(x0,A,B,R,QG,del2,M,G0.parametertype,opts);
    %x=udao_hessmult(x0,A,B,R,QG,del2,M,G0.parametertype,opts);
elseif strcmp(opt_method,'pso')
    x=udao_pso(x0,A,B,R,QG,del2,M,G0.parametertype,opts);
    x=udao_hessmult(x,A,B,R,QG,del2,M,G0.parametertype,opts);
elseif strcmp(opt_method,'nohess')
    x=udao_nohessian(x0,A,B,R,QG,del2,M,G0.parametertype,opts);
else
    x=udao_hessmult(x0,A,B,R,QG,del2,M,G0.parametertype,opts);
end

W.C=reshape(x(1:n*m),m, n);
W.D=reshape(x(m*n+1:m*n+m*m),m,m);
G=pfe(W,M);

if optget(opts,'plot_step3',0)
    plotXF(G,freq,Yv,[],'UDAO');
    if optget(opts,'plot_haeig',0)
        figure
        plot_haeig(G,[min(freq) max(freq)*10 -1 1],(max(freq)-min(freq))/1000);
    end
end











