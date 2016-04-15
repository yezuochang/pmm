function [err,L1,L2]=ss_error(S,F,H1,wgt_scheme)

if nargin < 4
    wgt_scheme = 5;
end

H2=ss_xf(S,F);
[err,L1,L2]=xf_error(H1,H2,wgt_scheme);

