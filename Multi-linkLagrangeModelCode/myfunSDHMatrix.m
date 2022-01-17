function Ai=myfunSDHMatrix(a,alpha,q,d)
%%
%Transform matrix in SDH coordinate 
%Ai,is the transform matrix between two coordinate
Ai=[cos(q),-sin(q)*cos(alpha), sin(q)*sin(alpha),a*cos(q);
    sin(q), cos(q)*cos(alpha),-cos(q)*sin(alpha),a*sin(q);
         0,        sin(alpha),        cos(alpha),       d;
         0,                 0,                 0,       1
   ];
     
end
%https://blog.csdn.net/easy_R/article/details/87939293