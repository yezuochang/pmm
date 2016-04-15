function [G, W, done] = LC (G0, W0, F, H, opts)
[G1,W1]=ASYM(G0, W0, F, H, opts);
if opts.parametertype == 'Y'
    [G,done]=LC_H(G1,F,H,opts);
    %G=lc_new(G1);
else
    [G,done]=LC_S(G1,F,H,opts);
end
W=[];

function [S2, done] = LC_H (S1, F1, Yv1, opts)

count = 0;
done = 0;
S2 = S1;
while ~done && count <= 40
    count = count+1;
    if optget(opts,'verbose',0) == 1
        fprintf('%d th iteration: \n',count);
    end
    [S2,done] = lc_sub_H(S2,opts);    
end
if ~done 
    warning('Warning : No convergence within 40 iterations ! ');
else
    if optget(opts,'verbose',0) == 1
        disp('Passivity was successfully enforced.');
    end
end

function [new_mdl,done]=lc_sub_H(mdl,opts)
method = 0;
global mmdl;
global result_dir;
global cktname;
mmdl = mdl;
new_mdl = mdl;
A=mdl.A;
B=mdl.B;
C=mdl.C;
D=mdl.D;
[n,m]=size(B);
I=eye(n);
bw_min=optget(opts,'bw_min',1e-5);
mag_min=optget(opts,'mag_min',1e-5);

%[M,N]=EHP_H(mdl);
% M=HamiltonianHybrid(mdl);
% r=eig(M);
% 
% r=r(find(abs(r)~=Inf));
% r2=violation(r);
% f2=imag(r2)/2/pi;    % Need verification!!!!!!!!!!
% f2=f2(find(f2>=0));

r2=passivity_violation(mdl);
f2=r2/2/pi;

% fprintf('violations:\n');
% fprintf('%e\n',f2);
if (length(f2)==0)
    done=1;
    new_mdl=mdl;
    return;
else
    done=0;
end

[M,N,J]=EHP_H(mdl);
rv=sort(2*pi*f2);
nv=length(rv);
V=eigvec(M,N,j*rv);
% calculate the step size
for c=1:nv
    v=V(:,c);
    v1=v(1:n);
    v2=v(n+1:2*n);
    v3=v(2*n+1:2*n+m);
    slope(c)=v3'*v3/imag(v2'*v1);
    
    if(0==slope(c))
        error('Error: slope infinite.');
    end
    slope(c)=1/slope(c);     % necessary???????
end
% if (slope(1)>0)    
%     error('Error: need to fix the dc passivity first.');    
% end

if (slope(end)<0)    
    figure;
    plot_haeig(mdl,[0 20 -1 1],1e-3,'r');
    error('Error: need to fix the asymptotic passivity first at f=Inf.');    
end

% calculating violation groups: denoting index 1,2,3,...
index = slope;
count = 0;
current_n = 1;
for c = nv:-1:1
    count = count + sign(slope(c));
    if 0 == count
        index(c) = current_n;
        current_n = current_n + 1;
    else
        index(c) = current_n;
    end
end
if count~=0
rv = [zeros(count,1);rv];       %adding zeros
slope = [-1*ones(1,count),slope];
nv = nv + count;
index = [current_n*ones(1,count),index];
else
    current_n = current_n-1;
end

% cmap = lines(current_n);
% figure;hold on;
%     plot_haeig1(mdl,[0.0125,1],'-k',[0.01,10,-1,3]);
%     for c = 1:current_n
%         ix = find(index == c);
%         scatter(rv(ix),zeros(length(ix),1),...
%             50,ones(length(ix),1)*cmap(c,:),'*');
%     end


vio = [];  % 4-by-m : \omega_1;\omega_2;\omega_0;magnitude
%===================method = 0 =================================
if 0 == method
    n_interp = 50;
    for c = 1:current_n
        ix = find(index == c);
        if length(ix)<=1
            fprintf('Potential bug here. Check it later.\n');            
            continue;
        end
        vio_header = rv(ix(1));
        vio_tail = rv(ix(end));
        step_interp = (vio_tail-vio_header)/n_interp;
        temp = [];
        for r_interp = vio_header:step_interp:vio_tail
            s=i*r_interp;
            H=C*((s*I-A)\B)+D;
            G=H+H';
            min_lambda=min(eig(G));
            temp = [temp,[r_interp;min_lambda]];
        end
        if isempty(temp)
            temp=temp;
        end
        [min_C,IX] = min(temp(2,:));
        vio = [vio,[vio_header;vio_tail;temp(1,IX);-min_C]];
    end
%===================method = 1 =================================
elseif 1 == method
    k=0;
%     loc=[];
%     bandwidth=[];
%     lmag=[]; 
    opts1 = optimset('Algorithm','interior-point');
    % forward looking
    for c=1:nv
        if slope(c)<0
            for d=c+1:nv
                if (slope(d)>0)
                    break;
                end
            end
            [rm,fval]=fmincon(@min_lamda,rv(c)+1e-5,[],[],[],[],rv(c),rv(d),[],opts1);
            newmin=1;
            for dd=1:size(vio,2)
                if (approximate(vio(3,dd),rm))
                    vio(1,dd) = min(vio(1,dd),rv(c));
                    vio(2,dd) = max(vio(2,dd),rv(d));
                    newmin=0;
                    break;
                end
            end
            if newmin
                vio = [vio,[rv(c);rv(d);rm;-min_lamda(rm)]];
            end        
        end
    end

    % backward looking
    for c=nv:-1:1
        if slope(c)>0
            for d=c-1:-1:1
                if (slope(d)<0)
                    break;
                end
            end
            [rm,fval]=fmincon(@min_lamda,rv(c)-1e-5,[],[],[],[],rv(d),rv(c),[],opts1);
            newmin=1;
            for dd=1:size(vio,2)
                if (approximate(vio(3,dd),rm))
                    vio(1,dd) = min(vio(1,dd),rv(d));
                    vio(2,dd) = max(vio(2,dd),rv(c));
                    newmin=0;
                    break;
                end
            end
            if newmin
                vio = [vio,[rv(d);rv(c);rm;-min_lamda(rm)]];
            end
        end
    end
%===================method = 2 =================================
else
    n_interp = 50;
    r_pair = [];   % 2-by-m matrix
    v_pair = [];   % ports-by-2-by-m 3D matrix
    for c = 1:current_n
        ix = find(index == c);
        if 2 == length(ix)
            vio_header = rv(ix(1));
            vio_tail = rv(ix(end));
            step_interp = (vio_tail-vio_header)/n_interp;
            temp = [];
            for r_interp = vio_header:step_interp:vio_tail
                s=i*r_interp;
                H=C*((s*I-A)\B)+D;
                G=H+H';
                min_lambda=min(eig(G));
                temp = [temp,[r_interp;min_lambda]];
            end
            [min_C,IX] = min(temp(2,:));
            vio = [vio,[vio_header;vio_tail;temp(1,IX);-min_C]];
        else
            ix1 = find(index == c & slope > 0);
            ix2 = find(index == c & slope < 0);
                if length(ix1)~=length(ix2)
                    warning('check point 002 : need debug');
                end
            temp_v2 = [];
            for d = 1:length(ix2)
                s=i*rv(ix2(d));
                H=C*((s*I-A)\B)+D;
                G=H+H';
                [VV,DD] = eig(G);
                [whatever,temp_ix] = min(abs(diag(DD)));
                temp_v2 = [temp_v2,VV(:,temp_ix)];
            end
            for d = 1:length(ix1)
                s=i*rv(ix1(d));
                H=C*((s*I-A)\B)+D;
                G=H+H';
                [VV,DD] = eig(G);
                [whatever,temp_ix] = min(abs(diag(DD)));
                temp_v1 = VV(:,temp_ix);
                [whatever,temp_ix] = max(abs(temp_v1'*temp_v2));
                r_pair = [r_pair,sort([rv(ix2(temp_ix));rv(ix1(d))])];
                v_pair = cat(3,v_pair,[temp_v2(:,temp_ix),temp_v1]);
                temp_v2(:,temp_ix) = 0*temp_v2(:,temp_ix);
            end
        end
    end
    % continuing calculating vio
    for c = 1:size(r_pair,2)
            vio_header = r_pair(1,c);
            vio_tail   = r_pair(2,c);
            step_interp = (vio_tail-vio_header)/n_interp;
            temp = [];
            for r_interp = vio_header:step_interp:vio_tail
                s=i*r_interp;
                H=C*((s*I-A)\B)+D;
                G=H+H';
                [VV,DD] = eig(G);
                temp1 = abs(v_pair(:,1,c)' * VV)+abs(v_pair(:,2,c)' * VV);
                [whatever,temp_ix] = max(temp1);
                min_lambda = DD(temp_ix,temp_ix);
                if min_lambda > 0
                min_lambda=min(eig(G));
                warning('check point 001 : need to debug');
                end
                temp = [temp,[r_interp;min_lambda]];
            end
            [min_C,IX] = min(temp(2,:));
            vio = [vio,[vio_header;vio_tail;temp(1,IX);-min_C]];
    end
end

%======================== plot  ================================
% if 0 < opts.plot
% n_vio = size(vio,2);
% cmap = lines(n_vio);
% figure;hold on;
%     plot_haeig1(mdl,[0.0125,1],'-k',[0.01,10,-1,3]);
%     for c = 1:n_vio
%         scatter(vio(1:3,c),[0;0;-vio(4,c)],...
%             50,ones(3,1)*cmap(c,:),'*');
%     end
% end
% if 2 == opts.plot
%     if ~exist(result_dir)
%     mkdir(result_dir);
%     end
%     get(0,'CurrentFigure');    
%     for try_num = 1:100
%         if ~exist(sprintf('%s/%s_eig%d.fig',result_dir,cktname,try_num))
%             break;
%         end
%     end
%     saveas(gcf,sprintf('%s/%s_eig%d.fig',...
%         result_dir,cktname,try_num));% potential problem with other compilers???
% end
%=================== compensate ================================
for c = 1:size(vio,2)
    if 0 == vio(2,c)
        a = 0;vio(3,c) = 0;
    elseif 0 == vio(1,c)
        a=1/sqrt(2)*abs(vio(3,c)^2/vio(2,c)-vio(2,c));
    else
        a=1/sqrt(2)*max(    abs( vio(3,c)^2/vio(1,c)-vio(1,c) ),...
            abs( vio(3,c)^2/vio(2,c)-vio(2,c) )    );
    end
    a = max(a, bw_min);
    b=vio(3,c)^2;
    K=max(1*a*vio(4,c),1e-3);
    I=eye(n);
    H=C*((i*vio(3,c)*I-A)\B)+D;
    G=H+H';
    [V1,D1]=eig(G);
    d1=diag(D1);
    [cc,ix]=min(d1);
    
    
  if(abs(b)>1e-8) 
      p=2;
    A2=[0,1;-b,-a];
    B2=zeros(2,m);
    B2(2,ix)=1;
    B2=B2*V1';
    C2=zeros(m,2);
    C2(ix,2)=K;
    C2=V1*C2;
    D2=zeros(m,m);
  else
      p=1;
      a=a/2;
      K=K/2;    
    A2=[-a];
    B2=zeros(1,m);
    B2(1,ix)=1;
    B2=B2*V1';
    C2=zeros(m,1);
    C2(ix,1)=K;
    C2=V1*C2;
    D2=zeros(m,m);
  end
    
% connect the two systems

    A=[A,zeros(n,p);zeros(p,n),A2];
    B=[B;B2];
    C=[C,C2];
    D=D+D2;
    n=n+p;
end


new_mdl.B=B;
new_mdl.A=A;
new_mdl.C=C;
new_mdl.D=D;

function y=approximate(x1,x2)
y=abs(x1-x2)<1e-3;

function [S2, done] = LC_S (S1, F1, Yv1, opts)

count = 0;
done = 0;
S2 = S1;
while ~done && count <= 40
    count = count+1;
    if optget(opts,'verbose',0) == 1
        fprintf('%d th iteration: \n',count);
    end
    [S2,done] = lc_sub_S(S2,opts);
end
if ~done 
    warning('Warning : No convergence within 40 iterations ! ');
else
    if optget(opts,'verbose',0) == 1
        disp('Passivity was successfully enforced.');
    end 
end


function [new_mdl,done]=lc_sub_S(mdl,opts)
global mmdl;
global result_dir;
global cktname;
%   if ~exist(opts)  %%%%%%%%%%%%%
%     opts.plot = 1;
%   end
mmdl = mdl;
new_mdl = mdl;
A=mdl.A;
B=mdl.B;
C=mdl.C;
D=mdl.D;
[n,m]=size(B);
I=eye(n);
bw_min=optget(opts,'bw_min',1e-5);
mag_min=optget(opts,'mag_min',1e-5);


% [M,N]=EHP_S(mdl);
% r=eig(M,N);
% r=r(find(abs(r)~=Inf));
% r2=violation(r);
% f2=imag(r2)/2/pi;
% f2=f2(find(f2>=0));
[r2,f2]=passivity_violation(mdl);
% fprintf('violations:\n');
% fprintf('%e\n',f2);
if (length(f2)==0)
    done=1;
    new_mdl=mdl;
    return;
else
    done=0;
end

rv=sort(2*pi*f2);
rv1 = [0;rv(1:length(rv)-1)];
mid_rv = (rv+rv1)/2;
I=eye(n);
for c = 1:length(mid_rv)
    H=C*((i*mid_rv(c)*I-A)\B)+D;
    sign_mid(c) = sign(1-max(svd(H)));
end
ix = find(sign_mid > 0);
index = ones(1,length(rv));
for c = length(ix):-1:1
    index(1:ix(c)-1) = index(1:ix(c)-1)+1;
end
if sign_mid(1) < 0
    rv = [0;rv];
    index = [index(1),index];
end


% cmap = lines(index(1));
% figure;hold on;
%     plot_svd(mdl,[0.001,10,0,2],2e-3,'-k');
%     for c = 1:index(1)
%         ix = find(index == c);
%         scatter(rv(ix),ones(length(ix),1),...
%             100,ones(length(ix),1)*cmap(c,:),'*');
%     end

    vio = [];  % 4-by-m : \omega_1;\omega_2;\omega_0;magnitude
    n_interp = 50;
    for c = 1:index(1)
        ix = find(index == c);
        vio_header = rv(ix(1));
        vio_tail = rv(ix(end));
        step_interp = (vio_tail-vio_header)/n_interp;
        temp = [];
        for r_interp = vio_header:step_interp:vio_tail
            s=i*r_interp;
            H=C*((s*I-A)\B)+D;
            temp = [temp,[r_interp;max(svd(H))]];
        end
        if isempty(temp)
            vio = [vio,[vio_header;vio_tail;mean(vio_header,vio_tail);1]];
            %%%%%%%%%%%%%%%%%%%%%%%%%   problem    %%%%%%%%%%%%%%%%%%%%%%%%
        else
            [max_C,IX] = max(temp(2,:));
            vio = [vio,[vio_header;vio_tail;temp(1,IX);max_C]];
        end
    end

%======================== plot  ================================
if 0 < opts.plot
n_vio = size(vio,2);
cmap = lines(n_vio);
figure;hold on;
    plot_svd(mdl,[0.001,10,0,2],2e-3,'-k');
    for c = 1:n_vio
        scatter(vio(1:3,c),[1;1;vio(4,c)],...
            50,ones(3,1)*cmap(c,:),'*');
    end
end
if 2 == opts.plot
    if ~exist(result_dir)
    mkdir(result_dir);
    end
    get(0,'CurrentFigure');    
    for try_num = 1:100
        if ~exist(sprintf('%s/%s_eig%d.fig',result_dir,cktname,try_num))
            break;
        end
    end
    saveas(gcf,sprintf('%s/%s_eig%d.fig',...
        result_dir,cktname,try_num));
end
%=================== compensate ================================
for c = 1:size(vio,2)
    if 0 == vio(2,c)
        a = 0;vio(3,c) = 0;
    elseif 0 == vio(1,c)
        a=1/sqrt(2)*abs(vio(3,c)^2/vio(2,c)-vio(2,c));
    else
        a=1/sqrt(2)*max(    abs( vio(3,c)^2/vio(1,c)-vio(1,c) ),...
            abs( vio(3,c)^2/vio(2,c)-vio(2,c) )    );
    end
    % a/2 is the bandwidth, we don't want it to be too small    
    a = max(a, bw_min);
    b=vio(3,c)^2;
%     K=max(0.6*a*vio(4,c),1e-3);
    beta = min(vio(4,c),1.05);
    beta = max(beta,1.00001);
    K=max(a - a/(beta*vio(4,c)),1e-3);
    I=eye(n);
    H=C*((i*vio(3,c)*I-A)\B)+D;
    G=H'*H;
    [V1,D1]=eig(G);
    d1=diag(D1);
    [cc,ix]=max(d1);
    
  if(abs(b)>1e-8) 
      p=2;
    A2=[0,1;-b,-a];
    B2=zeros(2,m);
    B2(2,ix)=1;
    B2=B2*V1';
    C2=zeros(m,2);
    C2(ix,2)=-K;
    C2=V1*C2;
    D2=eye(m);
  else
      p=1;
      a=a/2;      
      K=K/2;
    A2=[-a];
    B2=zeros(1,m);
    B2(1,ix)=1;
    B2=B2*V1';
    C2=zeros(m,1);
    C2(ix,1)=-K;
    C2=V1*C2;
    D2=eye(m);
  end
    
% connect the two systems : series connection

    A=[A,zeros(n,p);B2*C,A2];
    B=[B;B2*D];
    C=[D2*C,C2];
    D=D2*D;
    n=n+p;
end


new_mdl.B=B;
new_mdl.A=A;
new_mdl.C=C;
new_mdl.D=D;
new_mdl.parametertype='S';

