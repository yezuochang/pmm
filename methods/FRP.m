function [G, W, done] = FRP (G0,W0, F, H, opts)
[G1,W1]=ASYM(G0, W0, F, H, opts);
if opts.parametertype == 'Y'
    [G,done]=FRP_H(G1,F,opts);
else
    [G,done]=FRP_S(G1,F,opts);
end
W=[];

function [S2, done] = FRP_H (S1, F1, opts)

[S1.R,S1.poles] = ss2pr(S1.A,S1.B,S1.C);
%S1.E = S1.K;

%opts1.plot.s_pass=2*pi*i*linspace(0,10,1001).'; 
%opts1.plot.ylim=[-1 3]; 
opts1.Niter_out=50; 
opts1.Niter_in=5;   %%%%%%%  inner iteration
opts1.parametertype='Y';
opts1.outputlevel=0; %Min. output to screen
%opts1.plot=opts.plot;
[S2,bigSfit_passive,opts2]=RPdriver(S1,i*F1,opts1);

done = 1;  %????????????????????????????????

function [S2, done] = FRP_S (S1, F1, opts)


[S1.R,S1.poles] = ss2pr(S1.A,S1.B,S1.C);
%S1.E = S1.K;

%opts1.plot.s_pass=2*pi*i*linspace(0,10,1001).'; 
%opts1.plot.ylim=[0 2]; opts1.Niter_out=50; 
opts1.Niter_in=5;   %%%%%%%  inner iteration
opts1.parametertype='S';
opts1.outputlevel=0; %Min. output to screen
%opts1.plot=opts.plot;
[S2,bigSfit_passive,opts2]=RPdriver(S1,i*F1,opts1);


done = 1;  %????????????????????????????????











