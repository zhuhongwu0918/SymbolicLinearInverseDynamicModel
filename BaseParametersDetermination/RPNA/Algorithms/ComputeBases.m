function [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases(model, N, M) 
%% [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases(model, N, M)
% Computes system level bases for the parameter nullspace and its
% orthogonal complement. It also gives a minimal basis of unit vectors
% whose span is complementary to the parameter nullspace. Cell array inputs
% N and M are given as outputs of the RPNA algorithm.

    % Setup
    num_params = model.NB*10;
    perp_dim = 0; % dimension of the orthogonal complement to the parameter nullspace
    null_dim = 0; % dimension of the parameter nullspace
    perp_inds = cell(model.NB,1); % cell array containing column indicies of the perp_basis for each joint
    null_inds = cell(model.NB,1); % cell array containing column indicies of the null_basis for each joint
    
    for i = 1:model.NB
        dim_i = size(N{i},1);
     
        perp_inds{i} = perp_dim +1 : perp_dim + dim_i;
        null_inds{i} = null_dim +1 : null_dim + 10 - dim_i;
        
        perp_dim = perp_dim + dim_i;
        null_dim = null_dim + 10 - dim_i;
    end
    
    % Matrices whose colums will provide bases for each space
    Null_Basis    = zeros( num_params, null_dim );
    Perp_Basis    = zeros( num_params, perp_dim );
    Perp_Basis_sym    = sym(zeros( num_params, perp_dim ));
    Minimal_Basis = zeros( num_params, perp_dim );
    
    R    = cell(model.NB,1); % Holds the inerital transfer basis for each joint
    
    % Loop through one by at a time
    for i = 1:model.NB%第一个刚体，第i个刚体
        i_inds = parameter_inds(i); % Indicies of the inertial parameters for the body after joint i
        %整体机构的惯性参数的索引号
        parent = model.parent(i);   
        R{i} =  CleanMat( null(N{i}) );%R{i}是N{i}矩阵的零空间,N矩阵从RPNA algorithm算法来

        Minimal_Basis( i_inds , perp_inds{i} ) = M{i}; %M的元素 10*7
        
        % Undetectable inertial transfers from the child
        Null_Basis   ( i_inds , null_inds{i} ) = R{i};%Null_Basis（1:10，1:3）
%       size(R{i}) = 10     3,从R构造矩阵Null_Basis
        if parent > 0
            % X transforms velocities from parent to child
            X = model.Xtree{i};
            X_sym = model.Xtree_sym{i};
            
            % transforms parameters from child to parent
            %从子连杆到父连杆的坐标转换
            AX = Transform_Parameters(X);
            AX_sym = Transform_Parameters(X_sym);
            
            parent_inds = parameter_inds(parent);
            
            % -AX*R{i} gives the undetectable inertial transfers to the
            % parent
            %-AX*R{i}计算不可辨识的参数undetectable inertial transfers
            Null_Basis(parent_inds, null_inds{i} )= -AX*R{i};
%             R{i}是N{i}矩阵的零空间
            % See paper for recursive derivation of the basis for the
            % orthogonal complement of the parameter nullspace
            Perp_Basis(i_inds,:)                  = AX'*Perp_Basis(parent_inds,:);
            Perp_Basis_sym(i_inds,:)              = simplify(AX_sym'*Perp_Basis_sym(parent_inds,:));
        end
        
        % Set entries close to 0,1,-1 to 0,1,-1 so that the symbolic output
        % will be clean at the end.
        NN = CleanMat( rref(N{i})' );
        
        Perp_Basis(i_inds, perp_inds{i} )     = NN;
        Perp_Basis_sym(i_inds, perp_inds{i} ) = NN; 
    end
end

function inds = parameter_inds(i)
    inds = (1:10) + 10*(i-1);
end