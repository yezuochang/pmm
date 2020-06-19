function [F1,Yv1,scale]=datascale(F,Yv,opts)

m=size(Yv,1);
Yv1=permute(Yv,[3,1,2]);
s=j*2*pi*F;
for c=1:m
    for d=1:m
        %Yv1(:,c,d)=Yv1(:,c,d).*(1+s/1e8);
    end
end
Yv1=permute(Yv1,[2,3,1]);
if (min(F)==0)
    scale=max(F);%/100;
else
    %scale=sqrt(min(F)*max(F));
    scale=max(F);
end
F1=F/scale;
% if F1(1)==0 
%     if ~optget(opts,'enforceDC',0)
%         F1=F1(2:end);
%         Yv1=Yv1(:,:,2:end);
%     else
%         tol = optget(opts,'tol',0);
%         if 0 == tol
%             tol = 1e-6;
%         end
%         if opts.parametertype == 'Y'
%             [V,eigYv1] = eig(Yv(:,:,1)+Yv(:,:,1)');
%             eigYv1 = diag(eigYv1);
%             if min(eigYv1) < 1.2*tol
%                 ix = find(eigYv1 < 1.2*tol);
%                 eigYv1(ix) = 1.2*tol*ones(size(eigYv1(ix)));
%                 eigYv1 = diag(eigYv1);
%                 temp = V*eigYv1*V';
%                 Yv1(:,:,1) = (temp + Yv(:,:,1)-Yv(:,:,1)')/2;
%             end
%         else
%             [U, svdYv1,V] = svd(Yv(:,:,1));
%             svdYv1 = diag(svdYv1);
%             if max(svdYv1) > 1-1.2*tol
%                 ix = find(svdYv1 > 1-1.2*tol);
%                 svdYv1(ix) = 1-1.2*tol*ones(size(svdYv1(ix)));
%                 svdYv1 = diag(svdYv1);
%                 Yv1(:,:,1) = U*eigYv1*V';
%             end
%         end
%     end
% end