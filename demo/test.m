close all
addpath("..");
pmm_setup

opts=pmm_default;
opts.method="sdp";
opts.wgt_scheme = 4;
opts.enforceDC=1;
opts.fmax = 1e10;
opts.num_samples = 400;
opts.freqinterp = 'lin';
infiles = [
    "D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_doubler_spur_ind.s5p", 
    "D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_doubler_spur_ind_tx.s5p", 
    "D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_cko_core_emx.s8p",
    "D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_to_tx_line.s8p"];
%     "D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_doubler_topmetal.s44p"];
    
[G,W,F,H]=pmm("D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_cko_core_emx.s8p", 5, opts);
for i=1:size(infiles, 1)
    [G,W,F,H]=pmm(infiles(i), 8, opts);
    figure; plot_xf(G,F,H);legend("Data","Model","Error");
    figure; plot_haeig(G);
end
% infile="D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_doubler_topmetal.s44p";
% [G,W,F,H]=pmm(infile, 5, opts);
% figure; plot_xf(G,F,H);legend("Data","Model","Error");
% figure; plot_haeig(G);

% [G,W,F,H]=pmm("D:\workspace\pmm\data\huada\nportfiles\tc13\sp55_uniform.s64p", 2, opts);
% figure; plot_xf(G,F,H);legend("Data","Model","Error");
% figure; plot_haeig(G);