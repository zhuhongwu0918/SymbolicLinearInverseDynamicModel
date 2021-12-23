function A = thresholdsetzero(A,thres)%得到矩阵A的行号和列号
[m,n] = size(A);
for i=1:m
    for j=1:n
        if(A(i,j)<thres)
            A(i,j)=0;
        end
    end
end