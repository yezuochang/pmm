function [mdl1,P]=convertToComplexModel(mdl)

A=mdl.A;
B=mdl.B;
C=mdl.C;
D=mdl.D;

[P,L]=eig(full(A));

A=L;
B=P\B;
C=C*P;

% if (size(B,2)>1)
%     P1=diag(sum(B.').');
% else
%     P1=diag(B);
% end
% A=P1\(A*P1);
% B=P1\B;
% C=C*P1;

mdl1=mdl;
mdl1.A=A;
mdl1.B=B;
mdl1.C=C;



