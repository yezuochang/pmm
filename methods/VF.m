function [G,W]=VF(G0,W0,F,H,opts)
% VF -- generate state-space model from tabulated data using 
% Vector fitting method.
% This is simply an interface to VFdriver from Matrix Fitting Toolbox,
% This toolbox can be downloaded from 
% http://www.energy.sintef.no/produkt/VECTFIT/index.asp
%
% Copytright Zuochang Ye, 2012

opt1.N=opts.q ;%           %Order of approximation. 
opt1.Niter1=optget(opts,'vf_niter1', 200);    %Number of iterations for fitting sum of elements (fast!) --> Improved initial poles
opt1.Niter2=optget(opts,'vf_niter2', 100);    %Number of iterations for matrix fitting 
opt1.asymp=2;      %Fitting includes D   
opt1.logx=0;       %=0 --> Plotting is done using linear abscissa axis 
opt1.poletype=optget(opts,'poletype','lincmplx');
opt1.spy1=0; 
opt1.spy2=0; 
opt1.logx=0; 
opt1.logy=1; 
opt1.errplot= 1;  %%%%%%%%%%%%
opt1.phaseplot=0;
opt1.screen=0; %optget(opts,'verbose',0);
opt1.plot=opts.plot;
opt1.cmplx_ss=0;
opt1.weightparam=optget(opts,'wgt_scheme',3); %weight(s)=1/sqrt(abs(Hij(s)));
poles=[];
s1=2*pi*j*F;
Gc=VFdriver(H,s1,poles,opt1); 
Gc.parametertype=optget(opts,'parametertype','Y');

G=ss_real(Gc);

if F(1) == 0 && optget(opts,'enforceDC',0)
    [n,m] = size(G.B);
    Aeq = [kron(-G.B'*inv(G.A)',eye(m)),eye(m^2)];
    H0  = H(:,:,1);
    beq = vec(H0);
    F = F(2:end);
    H = H(:,:,2:end);
    [R,QG,del2]=preprocess(F,H,G,opts);    
    % y = R\QG;
    Q = R'*R;
    f = -QG'*R;    
    options.LargeScale = 'off';
    options.TypicalX = [vec(G.C);vec(G.D)];
    options.Display = 'off';
    y = quadprog(Q,f,[],[],Aeq,beq,[],[],[],options);
    y = real(y);
    G.C = reshape(y(1:m*n),m,n);
    G.D = reshape(y((m*n+1):(n*m+m^2)),m,m);
    G.D = H0 + G.C*(G.A\G.B);
end

W=[];
