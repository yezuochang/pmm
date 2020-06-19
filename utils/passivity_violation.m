function [r2,f2, slope]=passivity_violation(G,tol)
if nargin < 2
    tol = 0;
end
if tol < 0
    tol = 0;
end

% if G.parametertype == 'Y'
%     [r2,f2,slope]=half_size(G,tol);
% else
%     [r2,f2]=full_size(G,tol);
% end
[r2,f2]=full_size(G,tol);
f2=sort(f2);
r2=sort(r2);

function [r2,f2,slope]=full_size(G,tol);

if strcmp(G.parametertype,'Y')
    [M,N]=EHP_H(G);
else
    [M,N]=EHP_S(G);
end
r=eig(full(M),full(N));
ix=find(abs(real(r))*1e3<abs(imag(r)));
r1=imag(r(ix));
ix=find(r1>=0);
r2=r1(ix);
f2=r2/2/pi;

function [r2,f2,slope]=half_size(G,tol);
% 1. Adam Semlyen, A Half-Size Singularity Test Matrix for Fast and. Reliable Passivity Assessment
%    of Rational Models, 2009.
% 2. Bj?rn Gustavsen, Fast Passivity Assessment for S-Parameter Rational Models Via a 
%    Half-Size Test Matrix,2008

[A,B,C,D]=ss_data(G);
I=eye(size(D));
if strcmp(G.parametertype,'Y')
    M=A*(B*((D-0.5*tol*I)\C)-A);
    r=eig(M);
else
    gamma = 1-tol;
    M=(A-B*((D-gamma*I)\C))*(A-B*((D+gamma*I)\C));
    r=-eig(M);
end

r2=violation_halfsize(r);
[r2,slope]=violation_refine(G,r2,tol);
f2=r2/2/pi;

function [r2, slope]=violation_refine(G,r,tol)
[n,m] = size(G.B);
r2=[];
slope = [];

% for iter = 1:m
%     if G.parametertype == 'Y'
%         func=@(x) ss_lambda_min(G,x,iter);
%     else
%         func=@(x) ss_sigma_max(G,x,iter)-1;
%     end
%     options=optimset('TolFun',1e-16,'TolX',1e-6,'Display','off');
%     for c=1:length(r)
%         [x,fval,exitflag,output,jacobian]=fsolve(func,r(c),options);
%         if abs(fval) < 1e-6
% %         if exitflag >= 0
%             r2=[r2;x];
%             slope = [slope,jacobian];
%         end
%     end
% end
% 
% return;

% % fastest method
    if G.parametertype == 'Y'
        func=@(x) ss_lambda_near(G,x,tol);
    else
        func=@(x) ss_sigma_near(G,x,tol);
    end
    if 0 == tol
        tolfun = 1e-16;
    else
        tolfun = min(tol^2, 1e-16);
    end
    options=optimset('TolFun',tolfun,'TolX',1e-6,'Display','off');
    if 0 == tol
        F0 = 1e-6;
    else
        F0 = min(tol*0.1,1e-6);
    end
    for c=1:length(r)
        [x,fval,exitflag,output,jacobian]=fsolve(func,r(c),options);
        if abs(fval) < F0
%         if exitflag >= 0
            r2=[r2;x];
            slope = [slope,jacobian];
        end
    end
% merge same roots
if ~isempty(r2)
    [r2,ix] = sort(r2);
    slope = slope(ix);
    diff = r2 - [-1;r2(1:length(r2)-1)];
    diff = diff./r2;
    ix = find(diff > 1e-6);
    r2 = r2(ix);
    slope = slope(ix);
end
% find multiple roots
r1 = [];
slope1 = [];
for c = 1:length(r2)
    if G.parametertype == 'Y'
        [whatever, num] = ss_lambda_near(G,r2(c));
    else
        [whatever, num] = ss_sigma_near(G,r2(c));
    end
    r1 = [r1; repmat(r2(c),num,1)];
    slope1 = [slope1, repmat(slope(c),1,num)];
end
r2 = r1;
slope = slope1;

function r2=violation_halfsize(r)

%maxMag=max(abs(r));
%tol=maxMag*1e-6;
ix1=find(abs(real(r))>abs(imag(r)*1e8));
r1=r(ix1);
ix1=find(real(r1)>0);         %     why not equals zero ? 
r2=sort(sqrt(real(r1(ix1))));

