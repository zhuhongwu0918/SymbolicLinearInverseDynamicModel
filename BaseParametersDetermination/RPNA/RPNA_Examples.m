clear all; clc;
% addpaths

%% Options
% Uncomment to pick a model
% model = CreatePuma560();      % Industrial robot, 6 DoF
% model = CreateScara();        % Industrial robot, 3 DoF
model = CreateCheetahLeg();   % Single leg of a quadruped, 3 DoF
% model = CreateKirinLeg();        % Single leg of a quadruped, 3 DoF,RRP

MODEL_MOTORS = 0; % Include motor inertias (1), or ignore them (0)
FIXED_BASE   = 0; % Treat as fixed base (1), or ignore motion restrictions (0)

%model.gravity = [0 0 0]';
num_regressor_samples = 10; % For comparisson to numerical SVD

%% Compute Parameter Nullspace with SVD
fprintf(1,'\n\n***********************************************\n\n');
fprintf('Computing Random Regressors\n')
if MODEL_MOTORS
    Ystack = ComputeSampledRegressor_motor(model, num_regressor_samples);
    for i = 1:model.NB
        sub = Ystack(:,10*(model.NB+i-1)+1: 10*(model.NB+i-1)+10) ;
        Ystack(:,10*(model.NB+i-1)+1: 10*(model.NB+i-1)+10) = sub*model.motor_constraint{i};
    end
else
    Ystack = ComputeSampledRegressor(model, num_regressor_samples);
end
%随机采样生成观测矩阵Ystack

[Uy, Ey, Vy] = svd(Ystack);
%观测矩阵Ystack的SVD分解
if MODEL_MOTORS
    params_per_body = 20;
else
    params_per_body = 10;
end
SVD_Nullspace_Dimension = params_per_body*model.NB - rank(Ystack);

RangeBasis(1,1);
[~,SVD_Condition] = RangeBasis(Ystack');%RangeBasis给出一个与输入列跨度相同的满秩矩阵
% qr分解后正交矩阵Q的所有行和前满秩列的矩阵 重组会后的Ystack' [Q,R] =qr(Ystack'); 
%% Compute Parameter Nullspace with RPNA
%判断参数的可辨识性
RangeBasis(1,1);
param_names = {'m', 'mcx', 'mcy', 'mcz', 'Ixx', 'Iyy', 'Izz', 'Iyz', 'Ixz', 'Ixy'};
fprintf('Running RPNA\n');
fprintf(1,'\n===============================\n');
fprintf(1,'Identifiable Parameter Detail\n');

if MODEL_MOTORS 
    [N, M, V, C] = RPNA_motor(model,~FIXED_BASE);
    [~, RPNA_Condition] = RangeBasis(1);%%RangeBasis给出一个与输入列跨度相同的满秩矩阵
    [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases_motor(model, N, M);
    PrintParameterSummary_motor(model, N, M, V,C, param_names);
else%不考虑电机动力学，用以下形式
    [N, M, V, C] = RPNA(model,~FIXED_BASE);%N矩阵用来判定可辨识性
    [~, RPNA_Condition] = RangeBasis(1);
    [Null_Basis, Minimal_Basis, Perp_Basis, Perp_Basis_sym] = ComputeBases(model, N, M);
    %计算不可辨识参数的基，最小参数集的基等？？？
    %Perp_Basis_sym重要，用来求解惯性参数的线性组合系数，数值的Perp_Basis用来清洗矩阵
    PrintParameterSummary(model, N, M, V,C, param_names);
end

RPNA_Nullspace_Dimension = 0;
for i = 1:model.NB
    RPNA_Nullspace_Dimension = RPNA_Nullspace_Dimension + params_per_body-size(N{i},1);
end


%% Compute identifiable linear combinations with rref
%求线性组合和最小线性参数集合
fprintf(1,'===================================\n');
fprintf(1,'Minimal Parameter Detail\n');
fprintf(1,'===================================\n');
fprintf('Note: The listed linear cobminations of parameters are identifiable\n');
fprintf('from fully exciting data. These regroupings are also called minimal\n');
fprintf('parameters or base parameters in the literature.\n\n');

% Create variables for printing parameter regroupings
sym_params = sym( zeros(params_per_body*model.NB,1) ) ;    
%初始化标准惯性参数列向量，标准惯性参数*刚体个数
for i = 1:model.NB
    for k = 1:10
        sym_params(10*i-10 + k ) = sym(sprintf('%s%d',param_names{k},i));
        %生成标准参数集的符号表达式,m1 mcx1 mcy1 mcz1 Ixx1 Iyy1 Izz1 Iyz1 Ixz1 Ixy1 m2 mcx2...
        if MODEL_MOTORS
            sym_params(10*(i+model.NB)-10 + k ) = sym(sprintf('%sM%d',param_names{k},i));
        end
    end
end

% Compute identifiable parameter combinations from the basis for the
% subspace perpendicular to the parameter nullspace
Perp_Basis = rref(Perp_Basis')';
%自带函数rref，将矩阵转化为简化行阶梯形形式,Gauss-Jordan 消元法
Perp_Basis_sym = rref(Perp_Basis_sym')';
%将符号矩阵转化为简化行阶梯形形式,Gauss-Jordan 消元法

%靠近0的小值设为0
inds = find(abs(Perp_Basis) < 1e-8); % remove small values so printing is clean
Perp_Basis(inds) = 0;
Perp_Basis_sym(inds) = 0;
%靠近1的小值设为1
inds = find(abs(Perp_Basis-1) < 1e-8); % remove small values so printing is clean
Perp_Basis(inds) = 1;
Perp_Basis_sym(inds) = 1;
%靠近-1的小值设为-1
inds = find(abs(Perp_Basis+1) < 1e-8); % remove small values so printing is clean
Perp_Basis(inds) = -1;
Perp_Basis_sym(inds) = -1;

regrouping_matrix = sym(zeros(params_per_body*model.NB, params_per_body*model.NB  ));
%30*30的零矩阵，初始化
for i = 1:size(Perp_Basis_sym,2)%1:21,size(Perp_Basis_sym= 30    21
    ind = find(Perp_Basis_sym(:,i)==1,1);%Perp_Basis_sym第i列一列中等于1的,第一个索引
    regrouping_matrix(ind, :) = Perp_Basis_sym(:,i)';%将这一列换成行，放入regrouping_matrix
    
    % Identifable parameter combination 
    sym_result = Perp_Basis_sym(:,i)'*sym_params;%
    
    % Work to strip out zero coefficients去掉系数为0的项
    [coef, monomials] = coeffs(sym_result);
    %拆成两部分提取符号变量多项式的系数，monomials是符号变量表达式
    coef = CleanMat(coef);%设为0,1，-1，去掉过多的小数字位set values close to 0,1, and -1
    sym_result = simplify( sum( coef(:).*monomials(:) ) );
    
    % And then group terms that multiply each parameter
    sym_result = jacobian(sym_result, sym_params)*sym_params;
    %sym_result对标准参数求导，得到组合项，与2015年的J. Ros对应"Inertia transfer concept based general method for the determination of the base inertial parameters," Multibody System Dynamics, vol. 34, pp. 327-347, Aug .
    fprintf(1,'Regrouped parameter %s <= %s\n', char(sym_params(ind)), char(sym_result));
    %sym_result是线性组合后的关系
end

fprintf(1,'\n===================================\n');
fprintf(1,'Sanity Checks \n');
fprintf(1,'===================================\n');
fprintf('Null Check => norm( Ystack * Null_Basis ) = %e\n', norm( Ystack * Null_Basis , 'fro'))
fprintf('Perp Check => norm( Null_Basis\''*Perp_Basis ) = %e \n\n',norm(Null_Basis'*Perp_Basis,'fro') );

fprintf(1,'===================================\n');
fprintf(1,'Summary \n');
fprintf(1,'===================================\n');

if FIXED_BASE
    fprintf('Nullspace Dimension SVD  %d\n',SVD_Nullspace_Dimension)
end
fprintf('Nullspace Dimension RPNA %d\n',RPNA_Nullspace_Dimension)
fprintf('Identifiable Dimension %d\n',model.NB*params_per_body - RPNA_Nullspace_Dimension)

% fprintf('SVD Condition %f\n',SVD_Condition)
% fprintf('RPNA Condition %f\n',RPNA_Condition)

fprintf(1,'\n');