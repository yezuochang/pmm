function [G, W] = DAO (G0, W0, F, H, opts)
if opts.parametertype == 'Y'
    [G,W]=udao(G0,W0,F,H,opts);
else
    [G,W]=cdao(G0,W0,F,H,opts);
end