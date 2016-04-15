function G1 = enforce_dc(F,H,G,opts)

[R,QG,del2]=preprocess(F,H,G,opts);
[n,m] = size(G.B);
% y = R\QG;
x0 = vec(G.C);

Q = R'*R;
% Q = [Q1,Q2;Q2',Q3];
Q1 =  Q(1:m*n,1:m*n);
Q2 = Q(1:m*n,(m*n+1):(m*n+m^2));

f = vec(G.D)'*Q2'-QG'*R(:,1:m*n);
f = f';
Aeq = kron(-G.B'*inv(G.A)',eye(m));
beq = vec(H(:,:,1))-vec(G.D);

options.LargeScale = 'off';
options.TypicalX = x0;
options.Display = 'off';
y = quadprog(Q1,f,[],[],Aeq,beq,[],[],[],options);
% Q = R'*R;
% f = -QG'*R;
% Aeq = kron(-G.B'*inv(G.A)',eye(m));
% beq = vec(H(:,:,1))-vec(G.D);
% options.LargeScale = 'off';
% options.TypicalX = [vec(G.C);vec(G.D)];
% options.Display = 'off';
% y = quadprog(Q,f,[],[],Aeq,beq,[],[],[],options);
y = real(y);
G1 = G;
G1.C = reshape(y,m,n);
