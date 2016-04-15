function LegendreMATRIX=Legendre(X,N)
if N<2
    fprintf('N smaller than  needed parameter\n') ;%可能缺少一个跳出指令
end
    
LegendreMATRIX=ones(size(X,2),N);%生成Legendre多项式矩阵 
if N>=2
    LegendreMATRIX(:,2)=X';
end
if N>=3
    for i=3:N
    LegendreMATRIX(:,i)=[(2*(i-1)-1)*X(1,:)'.*LegendreMATRIX(:,i-1)-(i-2)*LegendreMATRIX(:,i-2)]/(i-1);%注意此时的阶次与矩阵序列号的关系，i与i-1的区别
    end
end