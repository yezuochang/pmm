function LegendreMATRIX=Legendre(X,N)
if N<2
    fprintf('N smaller than  needed parameter\n') ;%����ȱ��һ������ָ��
end
    
LegendreMATRIX=ones(size(X,2),N);%����Legendre����ʽ���� 
if N>=2
    LegendreMATRIX(:,2)=X';
end
if N>=3
    for i=3:N
    LegendreMATRIX(:,i)=[(2*(i-1)-1)*X(1,:)'.*LegendreMATRIX(:,i-1)-(i-2)*LegendreMATRIX(:,i-2)]/(i-1);%ע���ʱ�Ľ״���������кŵĹ�ϵ��i��i-1������
    end
end