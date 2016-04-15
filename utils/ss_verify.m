function ncss_verify(G,W,F,V,opts)
if optget(opts,'udao_verify',0) == 1
    udao_verify(G,W,F,V);
end