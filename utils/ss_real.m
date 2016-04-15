function [newMod,P]=ss_real(mod)
%% convert a model described by complex matrices to 
% a model described by pure real matrices.
% A is supposed to be a diagonal matrix.
% ref: Grivet, T-ADVP06
%%
if ss_isreal(mod)
    newMod=mod;
    P=[];
    return;
end
A=mod.A;
[n,n]=size(A);
dA=diag(A);
P=eye(n);
c=1;
while(c<=n)
    if (~isreal(dA(c)))
        P(c:c+1,c:c+1)=[1 -i;1 i]/sqrt(2);
        c=c+2;
    else
        c=c+1;
    end
end
newMod=mod;
newMod.A=real(P'*mod.A*P);
newMod.B=real(P'*mod.B);
newMod.C=real(mod.C*P);
newMod.D=real(mod.D);
if (isfield(mod,'K'))
    newMod.K=real(mod.K);
else
    newMod.K=real(mod.D*0);
end


