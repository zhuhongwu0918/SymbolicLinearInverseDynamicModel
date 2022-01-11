function bool = IsIdentifiable(model,N, i, k)
%% bool = IsIdentifiable(model,N, i, k)
% Determine if parameter k of body i is identifiable
% Requires the nullspace descriptor array N from the RPNA

    R = null(N{i});
    pi = zeros(10,1); %10个0，对应于10个惯性参数
    pi(k) = 1;%第k个设为1，如果k为1，则为m，k为2，则为mcx
    
    if norm( pi' * R ) > eps^.75
        bool = false;
        return
    end
    
    for j = 1:model.NB
        if i == model.parent(j)
            Rj = null(N{j});
            AX = Transform_Parameters(model.Xtree{j});
            %Xtree是6*6的旋转矩阵，表示两个坐标系间的变换关系
%             Transform_Parameters将6*6转换为10*10
            if norm( pi'*(AX*Rj) ) > eps^.75%R先乘一个坐标转换再判定
                bool=false;%如果这(AX*Rj)这一行有值,则不可辨识
                return
            end
        end
    end
    
    bool = true;
end