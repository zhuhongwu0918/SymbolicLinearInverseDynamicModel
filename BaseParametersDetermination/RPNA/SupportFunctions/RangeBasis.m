function [A, cA] = RangeBasis(A,clear)
%% Gives a full rank matrix with the same columnspan as the input

    persistent maxCond
    if isempty(maxCond)
        maxCond = 0;
    end
    if nargin == 2
        maxCond = 0;
    end
    [U E V] = svd(A);
    if min(size(E)) > 1
        E = diag(E);
    else
        if min(size(E)) > 0
            E = E(1,1);
        else
            E = [];
        end
    end
    r = sum( E > eps^.75 );
    if r == 0
        A = A(:,[]);
    else
        c = E(1) / E(r);%条件数，最大特征值除以最小特征值
        maxCond = max(c, maxCond);
        A = U(:,1:r);
    end
    cA = maxCond;

    [Q, ~,~] = qr(A);
    A = Q(:,1:rank(A));%Q是正交矩阵，R为上三角阵
end
    