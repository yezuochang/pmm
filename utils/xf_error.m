function [err,G1,G2]=xf_error(freq,H1,H2,opts)

wgt_scheme=optget(opts,'wgt_scheme',3);
if optget(opts,'freq_wgt_bw',0) == 0
    freq_weight = freq*0+1;    
else
    bw = optget(opts,'freq_wgt_bw',0);
    rho = freq(end)*bw;
    freq_weight = exp(-freq.*freq/(rho*rho));    
end
nf=size(H1,3);
m=size(H1,1);
m2=m*m;
err1=0;
for k=1:nf
    switch wgt_scheme
        case 1
            wgt=ones(m2,1);
        case 2
            wgt=1./vec(abs(H1(:,:,k)));            
        case 3
            wgt=1./vec(sqrt(abs(H1(:,:,k))));
        case 4
            wgt=ones(m2,1)./norm(H1(:,:,k));
        case 5
            wgt=ones(m2,1)./sqrt(norm(H1(:,:,k)));
        otherwise
            wgt=ones(m2,1);
    end
    wgt=freq_weight(k)*spdiags(wgt, 0, m2, m2);
    G1(:,k)=wgt*vec(real(H1(:,:,k)));
    G1(:,k+nf)=wgt*vec(imag(H1(:,:,k)));
    G2(:,k)=wgt*vec(real(H2(:,:,k)));
    G2(:,k+nf)=wgt*vec(imag(H2(:,:,k)));
%     (wgt*vec(imag(H1(:,:,k))))'
%     (wgt*vec(imag(H2(:,:,k))))'
    %err=err+vr'*vr+vi'*vi;
end

% err=norm(vec((G1-G2)))^2;
err=mean(vec(abs(G1-G2)));
