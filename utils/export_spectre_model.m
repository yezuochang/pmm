function export_spectre_model(G, W, subcktname, filename, mode)
if (~isreal(G.A))
    G=ss_real(G);
end
if (~isreal(W.A))
    W=ss_real(W);
end

G=ss_full(G);
W=ss_full(W);
[n,m]=size(G.B);

fid = fopen (filename,mode);
fprintf(fid,'subckt %s',subcktname);
fprintf(fid,' p%d n%d',[1:m;1:m]);
fprintf(fid,'\n');
fprintf(fid,'parameters A=[');
fprintf(fid,' %21.16e',reshape(G.A',n*n,1));
fprintf(fid,']\n');

fprintf(fid,'parameters B=[');
fprintf(fid,' %21.16e',reshape(G.B',n*m,1));
fprintf(fid,']\n');

fprintf(fid,'parameters C=[');
fprintf(fid,' %21.16e',reshape(G.C',n*m,1));
fprintf(fid,']\n');

fprintf(fid,'parameters D=[');
fprintf(fid,' %21.16e',reshape(G.D',m*m,1));
fprintf(fid,']\n');

fprintf(fid,'parameters Aw=[');
fprintf(fid,' %21.16e',reshape(W.A',n*n,1));
fprintf(fid,']\n');

fprintf(fid,'parameters Bw=[');
fprintf(fid,' %21.16e',reshape(W.B',n*m,1));
fprintf(fid,']\n');

fprintf(fid,'parameters Cw=[');
fprintf(fid,' %21.16e',reshape(W.C',n*m,1));
fprintf(fid,']\n');

fprintf(fid,'parameters Dw=[');
fprintf(fid,' %21.16e',reshape(W.D',m*m,1));
fprintf(fid,']\n');
fprintf(fid,'parameters noisetemp=27\n');
fprintf(fid,'parameters noiseref=2*P_K*(noisetemp+P_CELSIUS0)\n');
fprintf(fid,'\n');

fprintf(fid,'Xmod (');
for c=1:m
    if c==1
        fprintf(fid,'p%d n%d',c,c);
    else
        fprintf(fid,' p%d n%d',c,c);
    end
end
fprintf(fid,') cktrom a=A b=B c=C d=D\n');
for c=1:m
    fprintf(fid,'I%d (p%d n%d) cccs gain=-1 probe=V%d\n',c,c,c,c);
end

fprintf(fid,'Xnoise (');
for c=1:m
    if c==1
        fprintf(fid,'s%d 0',c);
    else
        fprintf(fid,' s%d 0',c);
    end
end
fprintf(fid,') cktrom a=Aw b=Bw c=Cw d=Dw\n');
for c=1:m
    fprintf(fid,'V%d (s%d 0) vsource  noisevec=[1 noiseref]\n',c,c);
end
fprintf(fid,'ends %s\n',subcktname);

fclose (fid);