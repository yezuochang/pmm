function [new_mdl,done]=EigPerturbForS(mod,alpha,K)% dC=(mod,r,alpha)

global mmdl;
mmdl=mod;
new_mdl=mod;
A=mod.A;
B=mod.B;
C=mod.C;
D=mod.D;
[n,m]=size(mod.B);
I = eye(n);

[M,N,J]=EHP_S(mod);
r=eig(M,N);
r=r(find(abs(r)~=Inf));
r2=violation(r);
f2=imag(r2)/2/pi;    % Need verification!!!!!!!!!!
f2=f2(find(f2>=0));
% fprintf('violations:\n');
% fprintf('%e\n',f2);
if (length(f2)==0)
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
                index1 = [index1,index(c)*ones(dif(c),1)];
            elseif dif(c) < 0
                rv1 = [rv1;rv(c)*ones(-dif(c),1)];
                slope1 = [slope1,+1e-8*ones(1,-dif(c))]; % small positive
                index1 = [index1,index(c)*ones(-dif(c),1)];
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

% calculating drv
drv = rv;
    n_interp = 50;
    for c = 1:current_n
        ix = find(index == c);
%         if 2 == length(ix)
%             vio_header = rv(ix(1));
%             vio_tail = rv(ix(end));
%             step_interp = (vio_tail-vio_header)/n_interp;
%             temp = [];
%             for r_interp = vio_header:step_interp:vio_tail
%                 s=i*r_interp;
%                 H=C*((s*I-A)\B)+D;
%                 G=H+H';
%                 min_lambda=min(eig(G));
%                 temp = [temp,[r_interp;min_lambda]];
%             end
%             [min_C,IX] = min(temp(2,:));
%             if vio_tail-vio_header > abs(min_C/slope(ix(1)))+abs(min_C/slope(ix(end)))
%                 drv(ix(1)) = abs(min_C/slope(ix(1)));
%                 drv(ix(end)) = -abs(min_C/slope(ix(end)));
%             else
%                 drv(ix(1)) = alpha*abs(vio_tail-vio_header);
%                 drv(ix(end)) = -alpha*abs(vio_tail-vio_header);
%             end
%         else
            for d = 1:length(ix)
                if slope(ix(d)) < 0 
                    drv(ix(d)) = alpha * (rv(ix(d-1))-rv(ix(d)));
                else
                    drv(ix(d)) = alpha * (rv(ix(d+1))-rv(ix(d)));
                end
            end
%         end
    end
%===================== plot =======================================
cmap = lines(current_n);
figure;hold on;
    plot_svd(mod,[0.001,15,0,2],2e-3,'-k');
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
if nargin < 3 || length(K) == 0
    K = eye(size(A));
end

%% solve the following quadratic programming 
% min ||dC K'||, s.t Zx=y   x = vec(dC)
%         Zk = Z (inv(K) kronecker I)
%         vec(dCk) = pinv(Zk)*y
%         dC = dCk*inv(K)'

Zk = Z*kron(inv(K),eye(length(D)));
xk = pinv(Zk)*y;      %      need verification
% xk = Zk\(drv*i);
dCk = reshape(xk,m,n);
dC = dCk * inv(K).';

% final perturbation
new_mdl.C = new_mdl.C+dC;


% display dC size
fprintf('Perturbation: |dCk| = %e\n',norm(dCk,'fro'));
fprintf('Perturbation: |dC| = %e\n',norm(dC,'fro'));

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



