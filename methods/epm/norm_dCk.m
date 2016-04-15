function f=norm_dCk(x)
global mm;
global nn;
dC=reshape(x,mm,nn);
global KK
% dCk=dC*KK';
dCk=dC;
f=norm(dCk,'fro'); %%%%%%%% dCk
