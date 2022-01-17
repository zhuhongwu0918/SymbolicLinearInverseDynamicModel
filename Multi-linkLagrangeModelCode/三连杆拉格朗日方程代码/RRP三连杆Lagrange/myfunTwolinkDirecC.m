function Cq=myfunTwolinkDirecC(n,D,q,dq)
%%仅针对二连杆适用
% syms q1 q2 
% syms dq1 dq2 
% q=[q1,q2]';
% dq=[dq1,dq2]';
% C=zeros(n,n);
% Ct=zeros(1,n)
%计算C(k,j)项
for k=1:n
    for j=1:n
        for i=1:n
%             Ct(i)=1/2*((diff(D(k,j),q(i))+(diff(D(k,i),q(j)))-(diff(D(i,j),q(k)))))*dq(i);
         Ct(i)=1/2*((diff(D(k,j),q(i))+(diff(D(k,i),q(j)))-(diff(D(i,j),q(k)))))*dq(i)*dq(j);
        end
        C(k,j)=sum(Ct);
    end
end
Cq=simplify(C);
 
 
end
%https://blog.csdn.net/easy_R/article/details/87939293