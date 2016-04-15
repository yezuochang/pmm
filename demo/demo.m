close all
addpath('..');
pmm_setup

opts=pmm_default;

%% A simple example to get start
% [G,W,F,H]=pmm('symind.s2p', 5, opts);
% plot_xf(G,F,H)

%% A medium-size problem

% Vector fitting only
opts.method='vf_only';
[G,W,F,H]=pmm('HHM1506.s3p', 20, opts);
figure; plot_xf(G,F,H);legend('Data','Model','Error');
figure; plot_haeig(G);

% SDP method 
opts.method='sdp';
[G,W,F,H]=pmm('HHM1506.s3p', 20, opts);
figure; plot_xf(G,F,H);

% Local compensation
opts.method='lc_only';
[G,W,F,H]=pmm('HHM1506.s3p', 20, opts);
figure; plot_xf(G,F,H);
%figure; plot_haeig(G);

% Local compensation + DAO method
opts.method='lc_dao';
[G,W,F,H]=pmm('HHM1506.s3p', 20, opts);
figure; plot_xf(G,F,H);legend('Data','Model','Error');
%figure; plot_haeig(G);