close all
addpath('..');
pmm_setup

opts=pmm_default;
opts.method='sdp';
opts.enforceDC=1;
% [G,W,F,H]=pmm('D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_doubler_spur_ind.s5p', 10, opts);
[G,W,F,H]=pmm('D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_cko_core_emx.s8p', 20, opts);
figure; plot_xf(G,F,H);legend('Data','Model','Error');
figure; plot_haeig(G);

% [G,W,F,H]=pmm('D:\workspace\pmm\data\huada\nportfiles\tc13\sp55_uniform.s64p', 2, opts);
% figure; plot_xf(G,F,H);legend('Data','Model','Error');
% figure; plot_haeig(G);   