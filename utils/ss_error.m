function [err,L1,L2]=ss_error(S,F,H1,opts)
H2=ss_xf(S,F);
[err,L1,L2]=xf_error(F,H1,H2,opts);

