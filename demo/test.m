close all
addpath('..');
pmm_setup

opts=pmm_default;
opts.method='sdp';
opts.enforceDC=1;
[G,W,F,H]=pmm('D:\workspace\pmm\data\huada\nportfiles\tc09\pll_sa_doubler_spur_ind.s5p', 10, opts);