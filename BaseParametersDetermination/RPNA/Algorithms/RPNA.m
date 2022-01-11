function [N, M, V , C] = RPNA(model, free_base)
%% [N, M, V , C] = RPNA(model, free_base)
% Recursive parameter nullspace algorithm
% The outputs of the algorithm are
% N: cell array for the nullspace descriptors characterizing undetectable
%    inertial transfers across each joint
% M: Basis for minimal paramers (sometimes called representative base
%    parameters) of each link
% V: Attainable velocity span basis for each link
% C: Null(C) gives the unidentifiable parameters of each link
%
% See the accompanying paper for further detail.
%


    if nargin == 1
        free_base = 0;
    end

    % Initialize Empty Cell arrays
    EmptyCells = cell(model.NB,1); 
    N   = EmptyCells;   C   = EmptyCells;   O   = EmptyCells;
    M   = EmptyCells;   V   = EmptyCells;   
    V_J = EmptyCells;  
    
    % Main algorithm loop
    for i =1:model.NB%从第i个刚体开始
        [~, Si] = jcalc( model.jtype{i}, 0 );%joint type，关节类型，Rz,Pz
        %返回旋转矩阵和关节自由度向量Si
        n_i = size(Si, 2);
        
        p = model.parent(i);
        if p == 0
            % If the base is to be treated as free, set the motion and
            % outer product bases to represent the full space. 
            % Else initialize them from the fixed base.
            if free_base
                Vp = eye(6);
                Cp = eye(10);
            else
                Vp =  get_gravity(model);
                Cp = zeros(0,10);
            end
        else
            Vp = V{p};
            Cp = C{p};
        end
        
        % Propage the velocity and outer product spans across the link
        V_J{i} = model.Xtree{i} * Vp;
       
        % Compute the collection of rate matricies for multi-DoF joints
        crmSet = cell( n_i, 1);
        paramRateSet     = cell( n_i, 1);
        for k= 1:size(Si,2)
            % Velocity rates of change with joint angle
            crmSet{k} = crm(Si(:,k) );%叉乘操作，6*6
            
            % Inertail parameters rates of change
            paramRateSet{k} = Rate_Parameters( Si(:,k) );
            %将自由度向量进行转换Rate_Parameters？？？原文中有描述
        end
        
        % Propagate the velocity across a joint  与算法中第4行对应
        V{i} = RangeBasis([ SwitchedControllabilityMatrix( crmSet , V_J{i} ) Si]); 
        %RangeBasis给出一个与输入列跨度相同的满秩矩阵
        O{i} = SwitchedObservabilityMatrix(Cp*Transform_Parameters( model.Xtree{i} ), paramRateSet );
        N{i} = OutputMatrix(V{i},Si);
        for k = 1:n_i
            N{i} = [N{i} ; O{i} * paramRateSet{k}];
        end
        N{i} = RangeBasis(N{i}')';
        C{i} = RangeBasis( [ Cp*Transform_Parameters( model.Xtree{i} ) ; N{i} ]')'; 
        M{i} = UnitVectorComplementarySubspace( null(N{i}) );
    end
end