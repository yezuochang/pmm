function S2=ss_y2s(S1,Z0)
% ss_y2s -- converting an S-parameter system to a y-parameter system
% Input
%   S1: ss system for s-parameter
%   Z0: reference impedance, default 50 
% Output
%   S2: ss system for y-parameter 
% Reference:
% Carlos P. Coelho, Joel R. Phillips, L. Miguel Silveira, 
% Passive Constrained Rational Approximation Algorithm using Nevanlinna-Pick
% Interpolation, DATE 2002.
if nargin < 2
    Z0=50;
end
assert(S1.parametertype == 'S');
[A,B,C,D]=ss_data(S1);
[n,m]=size(B);
I=eye(m);
T=inv(I+D);
Y1=1/sqrt(Z0);
S2.A=A-B*T*C;
S2.B=B*T*Y1;
S2.C=-2*Y1*T*C;
S2.D=Y1*(I-D)*T*Y1;
S2.parametertype='Y';
