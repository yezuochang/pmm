function S2=ss_y2s(S1,Z0)
% ss_y2s -- converting an S-parameter system to a y-parameter system
% Input
%   S1: ss system for y-parameter
%   Z0: reference impedance, default 50 
% Output
%   S2: ss system for s-parameter 
% Reference:
% Carlos P. Coelho, Joel R. Phillips, L. Miguel Silveira, 
% Passive Constrained Rational Approximation Algorithm using Nevanlinna-Pick
% Interpolation, DATE 2002.
if nargin < 2
    Z0=50;
end
assert(S1.parametertype == 'Y');
[A,B,C,D]=ss_data(S1);
[n,m]=size(B);
I=eye(m);
Z0=I*Z0;
Z1=sqrt(Z0);
Y0=inv(Z0);
Y1=inv(Z1);
T=inv(I+D);
T1=inv(I+Z1*D*Z1);
S2.A=A-B*inv(D+Y0)*C;
S2.B=2*B*Z1*T1;
S2.C=-T1*Z1*C;
S2.D=(I-Z1*D*Z1)*T1;
S2.parametertype='S';
