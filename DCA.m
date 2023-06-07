%DCA 用chol分解逆解线性方程组
function [xj,j]=DCA(G,H,Kmax,epsilon)
[m,n]=size(G);
xj=ones(n,1)/(m*n*m*n);%Dinbach x0
j=0;
R=chol(G);
L=R \ (R' \ H);
while(j<Kmax)
%     gxj=xj'*G*xj;
%     hxj=xj'*H*xj;
    xj1=xj;
    xjj=xj1/(m*n*m*n);
    h=xjj'*H*xjj;
    cj=(xjj'*G*xjj)/h;
    xj=DCA1(xj1,L,cj,100,1e-2);
    if norm(xj-xj1)<epsilon
        break;
    end
    j=j+1;
end
end

