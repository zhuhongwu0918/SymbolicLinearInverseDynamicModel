function c=myfunChristoffel(n,D,i,j,k)
%n,freedom of coordinate
%D,initial Matrix
%该函数仅仅针对两连杆模型
syms q1 q2 q3
q=[q1,q2, q3]';
c(i,j,k)=1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)));
end
%https://blog.csdn.net/easy_R/article/details/87939293