function [G, W, done] = EPM (G0,W0, F, H, opts)
[G1,W1]=ASYM(G0, W0, F, H, opts);
if opts.parametertype == 'Y'
    [G,done]=EPM_H(G1,F,H,opts);
else
    [G,done]=EPM_S(G1,F,opts);
    r=passivity_violation(G);
    if ~isempty(r)
        G.C=G.C*0.999;
        G.D=G.D*0.999;        
    end      
end
W=[];

function [S2, done] = EPM_H (S1, F1,H1,opts)

alpha=optget(opts,'epm_alpha',0.1);

% calculating K, default: use weighting scheme
if 0
K = lyapchol(S1.A,S1.B);
else
A1 = S1.A;
B1 = S1.B;
C1 = S1.C;
D1 = S1.D;
[A2,B2,C2,D2]=cheby2(30,50,max(F1)*3*pi,'low','s');
A2 = kron(eye(length(D1)),A2);
B2 = kron(eye(length(D1)),B2);
C2 = kron(eye(length(D1)),C2);
D2 = kron(eye(length(D1)),D2);
K = lyapchol([A1,B1*C2;zeros(length(A2),length(A1)),A2],[B1*D2;B2]);
K = K(1:length(A1),1:length(A1));
end
invK = inv(K);

% if optget(opts,'enforce_DC',0)
%     if F1(1) == 0
%         [R,QG,del2]=preprocess(F1(2:end),H1(:,:,2:end),S1,opts);
%     else
%         [R,QG,del2]=preprocess(F1,H1,S1,opts);
%     end
% else 
%     R = [];
% end

%%
count = 0;
done = 0;
S2 = S1;
tol = optget(opts,'tol',0);
if 0 == tol
    tol = 1e-6;
end
S2.D = S2.D - 0.55*tol*eye(length(S2.D));
while ~done && count <= 80
    count = count+1;
    if optget(opts,'verbose',0) == 1
        fprintf('%d th iteration: \n',count);
    end
    [S2,done] = eig_perturb_h(S2,alpha,invK,F1, H1,opts);
end
S2.D = S2.D + 0.55*tol*eye(length(S2.D));
if ~done 
    warning('Warning : No convergence within 40 iterations ! ');
else
    if optget(opts,'verbose',0) == 1
        disp('Passivity was successfully enforced.');
    end    
end

function [S2, done] = EPM_S (S1, F1, opts)

alpha=optget(opts,'epm_alpha',0.1);

% calculating K, default: use weighting scheme
if 0
K = lyapchol(S1.A,S1.B);
else
A1 = S1.A;
B1 = S1.B;
C1 = S1.C;
D1 = S1.D;
[A2,B2,C2,D2]=cheby2(30,50,max(F1)*3*pi,'low','s');
A2 = kron(eye(length(D1)),A2);
B2 = kron(eye(length(D1)),B2);
C2 = kron(eye(length(D1)),C2);
D2 = kron(eye(length(D1)),D2);
K = lyapchol([A1,B1*C2;zeros(length(A2),length(A1)),A2],[B1*D2;B2]);
K = K(1:length(A1),1:length(A1));
end
invK = inv(K);
%%
count = 0;
done = 0;
S2 = S1;
while ~done && count <= 40
    count = count+1;
    if optget(opts,'verbose',0) == 1
        fprintf('%d th iteration: \n',count);
    end
    [S2,done] = eig_perturb_s(S2,alpha,invK,opts);
end
if ~done 
    warning('Warning : No convergence within 40 iterations ! ');
else
    if optget(opts,'verbose',0) == 1
        disp('Passivity was successfully enforced.');
    end
end


function [new_mdl,done]=eig_perturb_h(mod,alpha,invK,F1,H1,opts)
new_mdl=mod;
A=mod.A;
B=mod.B;
C=mod.C;
D=mod.D;
[n,m]=size(mod.B);
I = eye(n);

[M,N,J]=EHP_H(mod);
[r2,f2,slope]=passivity_violation(mod);

if isempty(r2)
    done=1;
    new_mdl=mod;
    return;
else
    done=0;
end

rv=sort(2*pi*f2);

nv=length(rv);
V=eigvec(M,N,j*rv);

% test_V = rv';
% for c = 1:size(V,2)
%     test_V(c) = norm((M-i*rv(c)*N)*V(:,c))/norm(V(:,c));
% end
% test_V

% calculate the step size
% for c=1:length(rv)
%     v=V(:,c);
%     v1=v(1:n);
%     v2=v(n+1:2*n);
%     v3=v(2*n+1:2*n+m);
%     slope(c)=v3'*v3/imag(v2'*v1);
%     if(0==slope(c))
%         error('Error: slope infinite.');
%     end
%     slope(c)=1/slope(c);
%     flag(c)=sign(slope(c));
%     if 0 == flag(c)    %%%%%%%
%         disp('Warning : slope == 0 exists');
%         flag(c) = [];
%         slope(c) = [];
%         rv(c) = [];
%         V(:,c) = [];
%     end
% end

% slope = slope/2;      %???????????  necessary ? 
flag = sign(slope);
nv=length(rv);
% flag denotes the slope 
if (flag(end)== -1 )    
%     error('Error: need to fix the asymptotic passivity first.');    
    wrongslope = 1;
else
    wrongslope = 0;
end

% calculating violation groups: denoting index 1,2,3,...
index = slope;
count = 0;
current_n = 1;
if ~wrongslope
for c = nv:-1:1
    count = count + sign(slope(c));
    if 0 == count
        index(c) = current_n;
        current_n = current_n + 1;
    elseif count < 0
        wrongslope = 1;
    else
        index(c) = current_n;
    end
end
if count~=0
            temp = length(find(index == index(1))); %adding symetric roots
            rv    = [-rv(temp:-1:1);rv];
            index = [index((temp:-1:1)),index];
            slope = [-slope((temp:-1:1)),slope];
            V     = [V(:,(temp:-1:1)),V];
else
    current_n = current_n-1;
end
end % slope not wrong


% wrong slope: recalculating...
    if wrongslope
%         fprintf('Recalculating slopes... \n');'
        rv1 = [0;rv(1:length(rv)-1)];
        mid_rv = (rv+rv1)/2;
        I=eye(n);
        for c = 1:length(mid_rv)
            H=C*((i*mid_rv(c)*I-A)\B)+D;
            n_mid(c) = length(find(eig(H+H')<0));
        end
        % form index
        ix = find(n_mid == 0);
        index = ones(1,length(rv));
        for c = length(ix):-1:1
            index(1:ix(c)-1) = index(1:ix(c)-1)+1;
        end
        current_n = index(1);

        % calculating slope: positive or negative
        n_mid1 = [n_mid(2:end),0];
        dif = n_mid - n_mid1;
        
        slope1 = [];
        index1 = [];
        V1     = [];
        rv1    = [];
        for c = 1:length(dif)
            % remember rv, slope, index, V  AT THE SAME TIME!
            if dif(c) > 0
                rv1 = [rv1;rv(c)*ones(dif(c),1)];
                slope1 = [slope1,+1e-16*ones(1,dif(c))]; % small positive
                index1 = [index1,index(c)*ones(1,dif(c))];
                V1 = [V1,V(:,c)*ones(1,dif(c))];
            elseif dif(c) < 0
                rv1 = [rv1;rv(c)*ones(-dif(c),1)];
                slope1 = [slope1,-1e-16*ones(1,-dif(c))]; % small positive
                index1 = [index1,index(c)*ones(1,-dif(c))];
                V1 = [V1,V(:,c)*ones(1,-dif(c))];
            end
        end
        slope = slope1;
        index = index1;
        V     = V1;
        rv    = rv1;
        % adding symetric roots
        if n_mid(1) > 0
            temp = length(find(index == index(1)));
            rv    = [-rv(temp:-1:1);rv];
            index = [index((temp:-1:1)),index];
            slope = [-slope((temp:-1:1)),slope];
            V     = [V(:,(temp:-1:1)),V];
        end
        if isempty(rv)
             done=1;
            new_mdl=mod;
        return;
        else
            done=0;
        end
    end

V=eigvec(M,N,j*rv);   % only the first No.temp columns need this 
% calculating drv
drv = rv;
    n_interp = 50;
    for c = 1:current_n
        ix = find(index == c);
        if 2 == length(ix)
            vio_header = rv(ix(1));
            vio_tail = rv(ix(end));
            step_interp = (vio_tail-vio_header)/n_interp;
            temp = [];
            for r_interp = vio_header:step_interp:vio_tail
                temp = [temp,[r_interp;ss_lambda_min(mod,r_interp)]];
            end
            if ~isempty(temp)
                [min_C,IX] = min(temp(2,:));
                vio_mid = temp(1,IX);
            else  %  temp = [];
                vio_mid = (vio_tail+vio_header)/2;
                min_C = ss_lambda_min(mod,vio_mid);
            end
%             if vio_tail-vio_header > abs(min_C/slope(ix(1)))+abs(min_C/slope(ix(end)))
%                 drv(ix(1)) = abs(min_C/slope(ix(1)));
%                 drv(ix(end)) = -abs(min_C/slope(ix(end)));
%             else
%                 drv(ix(1)) = alpha*abs(vio_tail-vio_header);
%                 drv(ix(end)) = -alpha*abs(vio_tail-vio_header);
%             end
            if vio_mid-vio_header > abs(min_C/slope(ix(1)))
                drv(ix(1)) = abs(min_C/slope(ix(1)));
            else
                drv(ix(1)) = 2*alpha*abs(vio_mid-vio_header);
            end
            if vio_tail-vio_mid > abs(min_C/slope(ix(end)))
                drv(ix(end)) = -abs(min_C/slope(ix(end)));
            else
                drv(ix(end)) = -2*alpha*abs(vio_tail-vio_mid);
            end
        else
            for d = 1:length(ix)
                if slope(ix(d)) > 0 
                    drv(ix(d)) = alpha * (rv(ix(d-1))-rv(ix(d)));
                else
                    drv(ix(d)) = alpha * (rv(ix(d+1))-rv(ix(d)));
                end
            end
        end
    end
    
%========================= plot =================================
if optget(opts,'plot',0)
    cmap = lines(current_n);
    figure;hold on;
%     plot_haeig1(mod,[0.0125,1],'-k',[0.01,20,-1,3]);
  plot_haeig(mod);
    for c = 1:current_n
        ix = find(index == c);
        scatter(rv(ix)/2/pi,zeros(length(ix),1),...
            50,ones(length(ix),1)*cmap(c,:),'*');
        d = find(drv(ix)>0);
        scatter(rv(ix(d))/2/pi+drv(ix(d))/2/pi,zeros(length(ix(d)),1),...
            50,ones(length(ix(d)),1)*cmap(c,:),'>');
        d = find(drv(ix)<0);
        scatter(rv(ix(d))/2/pi+drv(ix(d))/2/pi,zeros(length(ix(d)),1),...
            50,ones(length(ix(d)),1)*cmap(c,:),'<');
    end
end
%================================================================
nv=length(rv);
% calculate the perturbation matrix
for c=1:nv
    v=V(:,c);
    v1=v(1:n);
    v3=v(2*n+1:2*n+m);
    y(c,1)=v'*J*N*v*drv(c)*i;
%     Z(c,:)=2*(kron(real(v1.'),real(v3.'))+...
%              kron(imag(v1.'),imag(v3.')));
      Z(c,:)=2*(kron(real(v1'),real(v3'))+...
               kron(imag(v1'),imag(v3')));
end

% set default K = eye
if nargin < 3 || isempty(invK)
    invK = eye(size(A));
end

% rv < 0 results in degenerated constrained problems
ix = find(rv<0);
y(ix) = [];
Z(ix,:) = [];

% add DC enforcement
if F1(1)==0 && optget(opts,'enforceDC',0)
    Aeq = kron(-B'*inv(A)',eye(m));
    beq = zeros(m^2,1);
    Z = [Z;Aeq];
    y = [y;beq];
end

% solve the following quadratic programming 
% min ||dC K'||, s.t Zx=d   x = vec(dC)
%         Zk = Z (inv(K) kronecker I)
%         vec(dCk) = pinv(Zk)*y
%         dC = dCk*inv(K)'

Zk = Z*kron(invK,eye(length(D)));
xk = pinv(Zk)*y;      %      need verification
% xk = Zk\(drv*i);
dCk = reshape(xk,m,n);
dC = dCk * invK.';

% final perturbation
new_mdl.C = new_mdl.C+dC;


% display dC size

if optget(opts,'verbose',0) == 1
    fprintf('Perturbation: |dCk| = %e\n',norm(dCk,'fro'));
    fprintf('Perturbation: |dC| = %e\n',norm(dC,'fro'));
end





function [new_mdl,done]=eig_perturb_s(mod,alpha,invK,opts)% dC=(mod,r,alpha)

new_mdl=mod;
A=mod.A;
B=mod.B;
C=mod.C;
D=mod.D;
[n,m]=size(mod.B);
I = eye(n);

[M,N,J]=EHP_S(mod);
% r=eig(M,N);
% r=r(find(abs(r)~=Inf));
% r2=violation(r);
% f2=imag(r2)/2/pi;    % Need verification!!!!!!!!!!
% f2=f2(find(f2>=0));
[r2,f2]=passivity_violation(mod);
% fprintf('violations:\n');
% fprintf('%e\n',f2);
if isempty(r2)
    done=1;
    new_mdl=mod;
    return;
else
    done=0;
end

rv=sort(2*pi*f2);

nv=length(rv);

        rv1 = [0;rv(1:length(rv)-1)];
        mid_rv = (rv+rv1)/2;
        I=eye(n);
        for c = 1:length(mid_rv)
            H=C*((i*mid_rv(c)*I-A)\B)+D;
            n_mid(c) = length(find(svd(H) > 1));
        end
        % form index
        ix = find(n_mid == 0);
        index = ones(1,length(rv));
        for c = length(ix):-1:1
            index(1:ix(c)-1) = index(1:ix(c)-1)+1;
        end
        current_n = index(1);

        % calculating slope: positive or negative
        n_mid1 = [n_mid(2:end),0];
        dif = n_mid - n_mid1;
        
        slope1 = [];
        index1 = [];
        rv1    = [];
        for c = 1:length(dif)
            % remember rv, slope, index, V  AT THE SAME TIME!
            if dif(c) > 0
                rv1 = [rv1;rv(c)*ones(dif(c),1)];
                slope1 = [slope1,-1e-8*ones(1,dif(c))]; % small negative
                index1 = [index1,index(c)*ones(1,dif(c))];
            elseif dif(c) < 0
                rv1 = [rv1;rv(c)*ones(-dif(c),1)];
                slope1 = [slope1,+1e-8*ones(1,-dif(c))]; % small positive
                index1 = [index1,index(c)*ones(1,-dif(c))];
            end
        end
        slope = slope1;
        index = index1;
        rv    = rv1;
        % adding symetric roots
        if n_mid(1) > 0
            temp = length(find(index == index(1)));
            rv    = [-rv(temp:-1:1);rv];
            index = [index((temp:-1:1)),index];
            slope = [-slope((temp:-1:1)),slope];
        end
        if (length(rv)==0)
             done=1;
            new_mdl=mod;
        return;
        else
            done=0;
        end
        
V=eigvec(M,N,j*rv);
nv = length(rv);
% calculating drv
drv = rv;
    n_interp = 50;
    for c = 1:current_n
        ix = find(index == c);
        if 2 == length(ix)
            vio_header = rv(ix(1));
            vio_tail = rv(ix(end));
            step_interp = (vio_tail-vio_header)/n_interp;
            temp = [];
            for r_interp = vio_header:step_interp:vio_tail
                temp = [temp,[r_interp;ss_sigma_max(mod,r_interp)]];
            end
            if ~isempty(temp)
                [max_C,IX] = max(temp(2,:));
                vio_mid = temp(1,IX);
            else  %  temp = [];
                vio_mid = (vio_tail+vio_header)/2;
                max_C = ss_sigma_max(mod,vio_mid);
            end
            drv(ix(1)) = 2*alpha*abs(vio_mid-vio_header);
            drv(ix(end)) = -2*alpha*abs(vio_tail-vio_mid);
        else
            for d = 1:length(ix)
                if slope(ix(d)) < 0 
                    drv(ix(d)) = alpha * (rv(ix(d-1))-rv(ix(d)));
                else
                    drv(ix(d)) = alpha * (rv(ix(d+1))-rv(ix(d)));
                end
            end
        end
    end
%===================== plot =======================================
if optget(opts,'plot',0)
    cmap = lines(current_n);
    figure;hold on;
    plot_svd(mod,[0.001,25,0,2],2e-3,'-k');
    for c = 1:current_n
        ix = find(index == c);
        scatter(rv(ix),ones(length(ix),1),...
            50,ones(length(ix),1)*cmap(c,:),'*');
        d = find(drv(ix)>0);
        scatter(rv(ix(d))+drv(ix(d)),ones(length(ix(d)),1),...
            50,ones(length(ix(d)),1)*cmap(c,:),'>');
        d = find(drv(ix)<0);
        scatter(rv(ix(d))+drv(ix(d)),ones(length(ix(d)),1),...
            50,ones(length(ix(d)),1)*cmap(c,:),'<');
    end
end
%===================================================================
    
% calculate the perturbation matrix
for c=1:nv
    v=V(:,c);
    v1=v(1:n);
    v4=v(2*n+m+1:2*n+2*m);
    y(c,1)=v'*J*N*v*drv(c)*i;
%     Z(c,:)=2*(kron(real(v1.'),real(v4.'))+...
%              kron(imag(v1.'),imag(v4.')));
      Z(c,:)=2*(kron(real(v1'),real(v4'))+...
               kron(imag(v1'),imag(v4')));
end

% set default K = eye
if nargin < 3 || isempty(invK)
    invK = eye(size(A));
end

%% solve the following quadratic programming 
% min ||dC K'||, s.t Zx=y   x = vec(dC)
%         Zk = Z (inv(K) kronecker I)
%         vec(dCk) = pinv(Zk)*y
%         dC = dCk*inv(K)'

Zk = Z*kron(invK,eye(length(D)));
xk = pinv(Zk)*y;      %      need verification
% xk = Zk\(drv*i);
dCk = reshape(xk,m,n);
dC = dCk * invK.';

% final perturbation
new_mdl.C = new_mdl.C+dC;


% display dC size

if optget(opts,'verbose',0) == 1
    fprintf('Perturbation: |dCk| = %e\n',norm(dCk,'fro'));
    fprintf('Perturbation: |dC| = %e\n',norm(dC,'fro'));
end

% 
% % to be deleted
% global mm;
% global nn;
% global KK;
% KK=K;
% mm=m;
% nn=n;
% x0=zeros(n*m,1);
% x=fmincon(@norm_dCk,x0,[],[],Z,y);
% fprintf('Perturbation: |dC| = %e\n',norm(dC,'fro'));
% 



