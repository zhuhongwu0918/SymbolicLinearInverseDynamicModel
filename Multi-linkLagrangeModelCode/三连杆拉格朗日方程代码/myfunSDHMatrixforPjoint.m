function Ai=myfunSDHMatrixforPjoint(a,alpha,q,d)
%%
%Transform matrix in SDH coordinate 
%Ai,is the transform matrix between two coordinate
Ai=[cospi(q),-sinpi(q)*cospi(alpha), sinpi(q)*sinpi(alpha),a*cospi(q);
    sinpi(q), cospi(q)*cospi(alpha),-cospi(q)*sinpi(alpha),a*sinpi(q);
         0,        sinpi(alpha),        cospi(alpha),       d;
         0,                 0,                 0,       1
   ];
     
end
%https://blog.csdn.net/easy_R/article/details/87939293