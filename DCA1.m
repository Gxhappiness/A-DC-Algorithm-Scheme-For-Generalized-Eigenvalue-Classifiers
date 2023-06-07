%DCA1 用chol分解逆解线性方程组
function xk=DCA1(xk,L,cj,Kmax,epsilon)
k=0;
while(k<Kmax)
    xk1=xk;
%     b=cj*xk;
%     xk=R \ (R' \ b);
    xk=L*cj*xk;
    if norm(xk-xk1)<epsilon
        break;
    end
    k=k+1;
end
end

