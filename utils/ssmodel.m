function S=ssmodel(A,B,C,D,parametertype)
if nargin < 5 
    type = 'H';
end
S.A = A;
S.B = B;
S.C = C;
S.D = D;
S.parametertype = parametertype;