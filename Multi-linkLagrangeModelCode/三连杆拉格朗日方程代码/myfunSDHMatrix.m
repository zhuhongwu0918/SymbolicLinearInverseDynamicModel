function Ai=myfunSDHMatrix(a,alpha,q,d)
%%
%Transform matrix in SDH coordinate 
%Ai,is the transform matrix between two coordinate
Ai=[cos(q),-sin(q)*cospi(alpha), sin(q)*sinpi(alpha),a*cos(q);
    sin(q), cos(q)*cospi(alpha),-cos(q)*sinpi(alpha),a*sin(q);
         0,        sinpi(alpha),        cospi(alpha),       d;
         0,                 0,                 0,       1
   ];
     
end
%https://blog.csdn.net/easy_R/article/details/87939293