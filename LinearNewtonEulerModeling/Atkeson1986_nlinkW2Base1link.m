clc;
clear;
close all

%计算基座以后的情况
%齐次变换矩阵from j to i,Z-Y-X(alpha, beta, gamma)
syms alpha beta gamma dalpha dbeta dgamma px0 py0 pz0  ddpx0 ddpy0 ddpz0 real
syms g real;
g_vec = [0 0 -g]';% 线加速度，以世界坐标为准
R_w_0 = [cos(alpha)*cos(beta),cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);...
    sin(alpha)*cos(beta),sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*cos(gamma);...
    -sin(beta),cos(beta)*sin(gamma),cos(beta)*cos(gamma)];
p_w_0 = [px0 py0 pz0]';
T_w_0 = [cell2sym([sym2cell(R_w_0),sym2cell(p_w_0)]);0, 0, 0,1];

z_w_0_w = eye(3,3);
% z_w_0_v = eye(3,3);
% z_w_0 = eye(6,6);
%从基座0开始迭代
w_w  = [0 0 0]';%角速度
dw_w = [0 0 0]';%角加速度
dv_w = g_vec;

syms ddalpha ddbeta  ddgamma real %角加速度和线加速度
w_0 = [dalpha dbeta dgamma]';
dw_0 = [ddalpha ddbeta ddgamma]';
w_0   = R_w_0*w_w + z_w_0_w*w_0;%公式1
dw_0  = R_w_0*dw_w + cross(R_w_0*w_w,z_w_0_w*w_0) + z_w_0_w*dw_0;
%仅使用于旋转关节
dv_0 = [ddpx0 ddpy0 ddpz0]' + dv_w;

%计算基座以后的情况
syms q1 q2 L1 L2 real
% modified DH = [alpha, a, d, theta]
dh_params = [0, 0, 0, q1;
             0, L1,0, 0;];%单连杆模型
[rows,~] = size(dh_params);
syms  m0 cx0 cy0 cz0 Ixx0  Ixy0  Ixz0  Iyy0  Iyz0  Izz0 real
Phibase = [m0,cx0*m0,cy0*m0,cz0*m0,  Ixx0,  Ixy0,  Ixz0,  Iyy0,  Iyz0,  Izz0]';
for i = 1:rows-1% 1     2  
    eval(['syms ','q',num2str(i),' real;']);%定义q1符号变量
    eval(['syms ','dq',num2str(i),' real;']);%定义dq1符号变量
    eval(['syms ','ddq',num2str(i),' real;']);%定义ddq1符号变量
    eval(['syms ','Lc',num2str(i),' real;']);%
    eval(['syms ','m',num2str(i),' real;']);%
    eval(['syms ','Ixx',num2str(i), ' Ixy',num2str(i), ' Ixz',num2str(i), ' Iyy',num2str(i), ' Iyz',num2str(i), ' Izz',num2str(i),' real;']);
    eval(['syms ','cx',num2str(i), ' cy',num2str(i), ' cz',num2str(i),' real;']);
    
    eval(['q(i)=','q',num2str(i),';']);%定义q(i)=q1符号向量
    eval(['dq(i)=','dq',num2str(i),';']);%定义dq(i)=dq1符号向量
    eval(['ddq(i)=','ddq',num2str(i),';']);%定义ddq(i)=ddq1符号向量
    eval(['mass_centerOi(i,:)=','[Lc',num2str(i),',0,0];']);%
    eval(['CoMOi(:,:,i)=','[cx',num2str(i),'; cy',num2str(i),'; cz',num2str(i),'];'])
    eval(['inertialOi(:,:,i)=','[Ixx',num2str(i), ' Ixy',num2str(i), ' Ixz',num2str(i), '; Ixy',num2str(i), ' Iyy',num2str(i), ' Iyz',num2str(i), '; Ixz',num2str(i), ' Iyz',num2str(i), ' Izz',num2str(i),'];']);%
    eval(['Phi(:,:,i)=','[m',num2str(i), ', m',num2str(i), '*cx',num2str(i), ', m',num2str(i), '*cy',num2str(i), ', m',num2str(i), '*cz',num2str(i), ',Ixx',num2str(i),', Ixy',num2str(i), ', Ixz',num2str(i), ', Iyy',num2str(i), ', Iyz',num2str(i), ', Izz',num2str(i) ']'';']);%
end

f_ext = [0 0 0 0 0 0]';%第一行是力；;%第二行是扭矩
%计算旋转矩阵和相邻关节位置
for j=1:rows-1
    %        z(:,j) = T(1:3,3,j);
    z(:,:,j) = [0 0 0 0 0 1]';%默认关节轴都在z方向上
end
z3 = [];%关节轴向的单位矢量，以关节坐标为参考系
%计算旋转矩阵和相邻关节位置
for j=1:rows-1
    %        z(:,j) = T(1:3,3,j);
    z3(:,:,j) = [0 0 1]';%默认关节轴都在z方向上
end
for i = 1:rows
    dh = dh_params(i,:);
    alpha(i) = dh(1);
    a(i) = dh(2);
    d(i) = dh(3);
    theta(i) = dh(4);
    if i == rows
        q(i) = 0;
    end
%   actuated axis transform matrix T(i-i -> i):
%MDH方法
   T(:,:,i) =  [cos(q(i))              ,-sin(q(i))              ,  0            ,  a(i);
                sin(q(i))*cos(alpha(i)), cos(q(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i);
                sin(q(i))*sin(alpha(i)), cos(q(i))*sin(alpha(i)),  cos(alpha(i)),  cos(alpha(i))*d(i);
                0                      , 0                      ,  0            ,  1                 ;];
            
   T_temp = T(:,:,i);
   R(:,:,i) = simplify(T_temp(1:3,1:3));
   P(:,:,i) = T_temp(1:3,4:4);%向量从i关节指向i+1关节
end

%基座初始值重新赋值
for i = 0:rows-1  %从基座0开始迭代
    if i == 0%定义基坐标，和重力方向
        wi  = w_0;%角速度
        dwi = dw_0;%角加速度
        dvOi = dv_0;
    else
        wi  = w(:,i);
        dwi = dw(:,i);
        dvOi= dvO(:,i);
    end
    if i+1 == rows
        break
    else
    w(:,:,i+1)   = R(:,:,i+1)*wi + dq(i+1)*z3(:,:,i+1);%公式1
    dw(:,:,i+1)  = R(:,:,i+1)*dwi + cross(R(:,:,i+1)*wi,dq(i+1)*z3(:,i+1)) + ddq(i+1)*z3(:,i+1);
    dvO(:,:,i+1)  = R(:,:,i+1)*( dvOi + cross(dwi,(p_w_0+P(:,:,i+1))) + cross(wi,cross(wi,(p_w_0+P(:,:,i+1)))));
    dvCi(:,:,i+1) = dvOi + cross(dwi,mass_centerOi(i+1,:)')...
        + cross(wi,cross(wi,mass_centerOi(i+1,:)'));
    end
end

for i = 1:rows-1
    f_i_i(:,:,i) = [dvO(:,:,i)-g_vec , crossMatrix(dw(:,:,i)) + crossMatrix(w(:,:,i))*crossMatrix(w(:,:,i)),zeros(3,6)]*Phi(:,:,i);
    n_i_i(:,:,i) = [zeros(3,1), crossMatrix((g_vec-dvO(:,:,i))), kMatrix(dw(:,:,i))+crossMatrix(w(:,:,i))*kMatrix(w(:,:,i))]*Phi(:,:,i);
    A_i(:,:,i) = [dvO(:,:,i)-g_vec , crossMatrix(dw(:,:,i)) + crossMatrix(w(:,:,i))*crossMatrix(w(:,:,i)),zeros(3,6);zeros(3,1), crossMatrix((g_vec-dvO(:,:,i))), kMatrix(dw(:,:,i))+crossMatrix(w(:,:,i))*kMatrix(w(:,:,i))];
    T_i(:,:,i) = [R(:,:,i)',zeros(3,3);crossMatrix(P(:,:,i))*R(:,:,i)',R(:,:,i)'];
    w_i_i(:,:,i) = [f_i_i(:,:,i);n_i_i(:,:,i)];
end
isequal(A_i(:,:,1)*Phi(:,:,1),w_i_i(:,:,1))

%the torque at joint 0 due to motion of base link （0） alone
f_0_0 = [dv_0(:,:,i)-g_vec , crossMatrix(dw_0) + crossMatrix(w_0)*crossMatrix(w_0),zeros(3,6)]*Phibase;
n_0_0 = [zeros(3,1), crossMatrix((g_vec-dv_0)), kMatrix(dw_0)+crossMatrix(w_0)*kMatrix(w_0)]*Phibase;
A_0 = [dv_0(:,:,i)-g_vec , crossMatrix(dw_0) + crossMatrix(w_0)*crossMatrix(w_0),zeros(3,6);zeros(3,1), crossMatrix((g_vec-dv_0)), kMatrix(dw_0)+crossMatrix(w_0)*kMatrix(w_0)];
T_0 = [R_w_0',zeros(3,3);crossMatrix(p_w_0)*R_w_0',R_w_0'];
w_0_0 = [f_0_0;n_0_0];

%构造连杆的完整U矩阵
for i = 1 : rows-1
    for j=rows-1:-1:i
        if j == i
            U(:,:,i,j) = A_i(:,:,i);
        else
            temp = 1;
            for k = i : j
                temp = temp * T_i(:,:,k);
            end
            U(:,:,i,j) = temp* A_i(:,:,i);
        end
    end
end
%%
%构造底座加单连杆的U矩阵
U_0_0 = A_0;
U_0_1 = T_0*T_i(:,:,1)*A_0;
%构造底座加单连杆的完整Y矩阵
w0 = [U_0_0,U_0_1]*[Phibase;Phi(:,:,1)];
w1 = [zeros(6,10),U(:,:,1,1)]*[Phibase;Phi(:,:,1)];
Y_one = [U_0_0,U_0_1;zeros(6,10),U(:,:,1,1)];
w = [w0;w1];
isequal(Y_one*[Phibase;Phi(:,:,1)],w)
% K_one = [z_w_0'*U_0_0,z'*U_0_1;z'*zeros(6,10),z'*U(:,:,1,1)];
%%
%按照测量方式向某自由度投影
for i = 1:size(U,3)
    for j = 1:size(U,4)
        K_tau(:,:,i,j) = z(:,:,i)' * U(:,:,i,j);
    end
end
K_tau_w_0 = [U_0_0,U_0_1];

K_tau_2d=cell(size(K_tau,3),size(K_tau,4));
K_tau_2d{1,1} = U_0_0;
K_tau_2d{1,2} = U_0_1;
K_tau_2d{2,1} = zeros(size(z(:,:,1)'*U_0_0));
K_tau_2d{2,2} = K_tau(:,:,1,1);
K_tau_2d = cell2sym(K_tau_2d);
%first 6 row is the base link, the last row is the first link.
%%
%赋值测试
ddpx0 = 2;ddpy0 = 0;ddpz0 = 0;
alpha=0; beta=0;   gamma=0; 
dalpha=0; dbeta=0;   dgamma=0; 
ddalpha=0; ddbeta=0;   ddgamma=0; 
q1=0.5*pi;dq1=0;ddq1=0;
px0 =1;py0=1;pz0=1;%基座距离世界坐标的位置
eval(K_tau_2d)
ans=thresholdsetzero(ans,1e-3);
ans*[Phibase;Phi(:,:,1)]