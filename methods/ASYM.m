function [G, W, done] = ASYM (G0, W0, F1, H1, opts)

if G0.parametertype == 'Y'
    [G]=ASYM_H(G0,F1,H1,opts);
else
    [G]=ASYM_S(G0,F1,H1,opts);
end
W=[];

function S1=ASYM_H(S1,F1,Yv1,opts)
[V,D] = eig((S1.D+S1.D')/2);
d = diag(D);
if min(d) <= 0
    if optget(opts,'verbose',0) == 1
        disp('Fixing asymptotical passivity... ');
    end
    
    ix = find(d < 0);
    d(ix) = -0.1*d(ix);
    ix = find(d == 0);
    d(ix) = 0.01*ones(size(d(ix)));
    D = diag(d);
    D = V*D*V';
    DD = repmat(D,[1,1,length(F1)]);
    Yv_temp = Yv1 - DD;
    [n,m]=size(S1.B);
    opts.N=round(n/m) ;%           %Order of approximation.
    opts.poletype='linlogcmplx'; %Mix of linearly spaced and logarithmically spaced poles
    opts.Niter1=20;    %Number of iterations for fitting sum of elements (fast!) --> Improved initial poles
    opts.Niter2=10;    %Number of iterations for matrix fitting
    opts.logx=0;       %=0 --> Plotting is done using linear abscissa axis
    opts.poletype='lincmplx';
    opts.spy1=0;
    opts.spy2=0;
    opts.logx=0;
    opts.logy=1;
    opts.errplot=0;
    opts.phaseplot=0;
    opts.screen=0; %optget(opts,'verbose',0);
    opts.asymp=1;      %Fitting D = 0
    s1=2*pi*j*F1;
    poles=[]; 
    Sc=VFdriver(Yv_temp,s1,poles,opts);
    Sc.D = D;
    Sc.parametertype=S1.parametertype;
    S1=ss_real(Sc);
    if F1(1) == 0 && optget(opts,'enforceDC',0)
        S1 = enforce_dc(F1,Yv1,S1,opts);
    end
end

function S1=ASYM_S(S1,F1,Yv1,opts)

[U,D,V] = svd(S1.D);
d = diag(D);
if max(d) >= 1
    if optget(opts,'verbose',0) == 1
        disp('Fixing asymptotical passivity... ');
    end
    ix = find(d > 1);
    d(ix) = 1./(d(ix)).^0.5;
    ix = find(d == 1);
    d(ix) = 0.9*ones(size(d(ix)));
    D = diag(d);
    D = U*D*V';
    DD = repmat(D,[1,1,length(F1)]);
    Yv_temp = Yv1 - DD;
    [n,m]=size(S1.B);
    opts.N=round(n/m) ;%           %Order of approximation.
    opts.poletype='linlogcmplx'; %Mix of linearly spaced and logarithmically spaced poles
    opts.Niter1=20;    %Number of iterations for fitting sum of elements (fast!) --> Improved initial poles
    opts.Niter2=10;    %Number of iterations for matrix fitting
    opts.logx=0;       %=0 --> Plotting is done using linear abscissa axis
    opts.poletype='lincmplx';
    opts.spy1=0;
    opts.spy2=0;
    opts.logx=0;
    opts.logy=1;
    opts.errplot=0;
    opts.phaseplot=0;
    opts.screen=optget(opts,'verbose',0);
    opts.asymp=1;      %Fitting D = 0
    s1=2*pi*j*F1;
    poles=[]; 
    
    Sc=VFdriver(Yv_temp,s1,poles,opts);
    Sc.D = D;
    S1=ss_real(Sc);
        if F1(1) == 0 && optget(opts,'enforceDC',0)
        S1 = enforce_dc(F1,Yv1,S1,opts);
        end
    S1.parametertype='S';
end