clc;
clear;

syms q1 q2 L1 L2 real
% modified DH = [alpha, a, d, theta]
% dh_params = [0, 0, 0, q1;
%              0, L1,0, 0;];%单连杆模型
dh_params = [0, 0, 0, q1;
             0, L1,0, q2;
             0, L2,0, 0];%双连杆模型

[rows,~] = size(dh_params);
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
z = [];%关节轴向的单位矢量，以关节坐标为参考系
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

syms g real;
g_vec = [g,0 0]';% 线加速度，以基坐标为准，
% 第一个关节会体现在旋转矩阵的dh参数的alpha上，基坐标x向下，则为g00。
% 第一根杆是上下摆动，而不是平面摆动。
for i = 0:rows-1  %从基座0开始迭代
    if i == 0%定义基坐标，和重力方向
        wi  = [0 0 0]';%角速度
        dwi = [0 0 0]';%角加速度
        dvOi = g_vec;
%       dvOi = [0,0 -g]';
    else
        wi  = w(:,i);
        dwi = dw(:,i);
        dvOi= dvO(:,i);
    end
    if i+1 == rows
        break
    else
    %此处用R’，为什么是反着的？
    w(:,:,i+1)   = R(:,:,i+1)*wi + dq(i+1)*z3(:,:,i+1);%公式1
    dw(:,:,i+1)  = R(:,:,i+1)*dwi + cross(R(:,:,i+1)*wi,dq(i+1)*z3(:,i+1)) + ddq(i+1)*z3(:,i+1);
    %dvO(:,:,1)=[g*cos(q1),g*sin(q1),0]^T
    dvO(:,:,i+1)  = R(:,:,i+1)*( dvOi + cross(dwi,P(:,:,i+1)) + cross(wi,cross(wi,P(:,:,i+1))));
    % dvO(:,:,i+1)  = R(:,:,i+1)*(dvOi + crossMatrix(dwi)*P(:,:,i+1) + cross(wi,cross(wi,P(:,:,i+1))));
    %序号是否有问题,--- ref : https://zhuanlan.zhihu.com/p/109770180
    %aci = dwi×rCi^i + wi×(wi×rCi^i) + aOi;
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

% for i = rows-1:-1:1
%     if i == rows-1
%         w(i) = w_i_i(:,:,i) + f_ext;
%     else
%         w(i) = w_i_i(:,:,1) + T_i(:,:,2)*w_i_i(:,:,2);
%     end
% end
%w_1_1 = A_1*Phi_1
%w_2_2 = A_2*Phi_2
isequal(A_i(:,:,1)*Phi(:,:,1),w_i_i(:,:,1))
% isequal(A_i(:,:,2)*Phi(:,:,2),w_i_i(:,:,2))
% T2 = [R(:,:,2),zeros(3,3);crossMatrix(P(:,:,2))*R(:,:,2),R(:,:,2)];
% w_2 = w_i_i(:,:,2) + f_ext;
% w_1_2 = T2*w_2_2;
% w_1 = w_1_1 + w_1_2;
% w_1 = w_i_i(:,:,1) + T_i(:,:,2)*w_i_i(:,:,2);

% tau_1 = z(:,:,1)'*w_1;
% tau_2 = z(:,:,2)'*w_2;

for i = 1 : rows-1
    for j=rows-1:-1:i
        if j == i
            U(:,:,i,j) = A_i(:,:,i);
        else
            temp = 1;
            for k = i : j
                temp = temp * T_i(:,:,k);
            end
%             fprintf('i=');i
%             fprintf('j=');j
            U(:,:,i,j) = temp* A_i(:,:,i);
        end
    end
end

for i = 1:size(U,3)
    for j = 1:size(U,4)
        K_tau(:,:,i,j) = z(:,:,i)' * U(:,:,i,j);
    end
end

K_tau_2d=cell(size(K_tau,3),size(K_tau,4));
for i = 1:size(K_tau,3)
     for j = 1:size(K_tau,4)
    K_tau_2d{i,j} =  K_tau(:,:,i,j);
     end
end
K_tau_2d = cell2sym(K_tau_2d);

% K = [U(:,:,1,1), U(:,:,1,2);zeros(6,10),U(:,:,2,2)];
% K_tau = [z(:,:,1)'*U(:,:,1,1), z(:,:,1)'*U(:,:,1,2);z(:,:,2)'*zeros(6,10),z(:,:,2)'*U(:,:,2,2)]
% 
% U_1_1 = A_i(:,:,1);
% U_1_2 = T_i(:,:,1)*T_i(:,:,2)*A_i(:,:,1);
% U_2_2 = A_i(:,:,2);
% % w_1 = [U_1_1,U_1_2]*[Phi_1;Phi_2];
% 
% K = [U_1_1, U_1_2;zeros(6,10),U_2_2];
% % K1 = [0,0,0,0,0,1] * K(0*6+1:1*6,:)
% % K2 = [0,0,0,0,0,1] * K(1*6+1:2*6,:)
% 
% K_tau = [z(:,:,1)'*U_1_1, z(:,:,1)'*U_1_2;z(:,:,2)'*zeros(6,10),z(:,:,2)'*U_2_2]
% % isequal(K_tau,[w_1;w_2])

%ID
load Vertical2link
g=9.81;
L1 = 0.2;
LC1 = 0.1;
L2 = 0.1;
LC2 = 0.05;
StartDataindex = 1000;
EndDataindex = 5000;
step = 10;
numOfData = (EndDataindex-StartDataindex)/step+1;

% K1_num=[];K2_num=[];K=[];tau=[];
for n=1:1:numOfData
    q1=getdatasamples(q1_ts,step*(n-1)+StartDataindex);
    dq1=getdatasamples(w1_ts,step*(n-1)+StartDataindex);
    ddq1=getdatasamples(b1_ts,step*(n-1)+StartDataindex);
    tau1_actual(n,1) = getdatasamples(t1_ts,step*(n-1)+StartDataindex);
%     tau(:,:,n) = [tau1_actual(n,1)];
    %二连杆
    q2=getdatasamples(q2_ts,step*(n-1)+StartDataindex);
    dq2=getdatasamples(w2_ts,step*(n-1)+StartDataindex);
    ddq2=getdatasamples(b2_ts,step*(n-1)+StartDataindex);
    tau2_actual(n,1) = getdatasamples(t2_ts,step*(n-1)+StartDataindex);
    tau_num(n,:) = [tau1_actual(n,1);tau2_actual(n,1)];
    tau(:,:,n) = [tau1_actual(n,1);tau2_actual(n,1)];
     
      K_num(:,:,n) = eval(K_tau_2d);
end

K_num_vec=[];tau_vec=[];
for i = 1:size(K_num,3)
    K_num_vec=[K_num_vec;K_num(:,:,i)];
    tau_vec = [tau_vec;tau(:,:,i)];
end
% K_num_vec =reshape(K_num,[401*2,20]);
% K_num_vec =cat(3,K_num);
% tau_vec=squeeze(K_num);

%自带pinv函数
Phi_cal = pinv(K_num_vec)*tau_vec
%%

Phi_cal = pinv([K1_num;K2_num])*[tau1_actual;tau2_actual]
Phi2_cal = pinv(K2_num)*tau2_actual
% [m1, m1*cx,m1*cy,m1*cz,Ixx1,Ixy1,Ixz1,Iyy1,Iyz1,Izz1]'
%          0
%     0.1032
%     0.0014
%          0
%          0
%          0
%          0
%          0
%          0
%     0.0144
%%
%圆柱体几何计算求解验证Ip
cx1=0.1;cy1=0;cz1=0;
cx2=0.05;cy2=0;cz2=0;
R_p_to_q1 = [cos(-0.5*pi),0,sin(-0.5*pi);0,1,0;-sin(-0.5*pi),0,cos(-0.5*pi)];
R_p_to_q2 = [cos(-0.5*pi),0,sin(-0.5*pi);0,1,0;-sin(-0.5*pi),0,cos(-0.5*pi)];
%
p1=[cx1;0;0];
p2=[cx2;0;0];
cx1_sim=0;cy1_sim=0;cz1_sim=0;
cx2_sim=0;cy2_sim=0;cz2_sim=0;
CoM1 = p1 + [cx1_sim;cy1_sim;cz1_sim];
CoM2 = p2 + [cx2_sim;cy2_sim;cz2_sim];
%
m1 =0.0628319;
m2 =0.0314159;
XX1=2.1090e-04;YY1=2.1090e-04;ZZ1= 3.1400e-06;
XX2=2.69653e-05;YY2=2.69653e-05;ZZ2= 1.5708e-06;
principalmoments1 = [XX1,0,0;0,YY1,0;0,0,ZZ1];
principalmoments2 = [XX2,0,0;0,YY2,0;0,0,ZZ2];
%
Iq1 = R_p_to_q1 * principalmoments1 * R_p_to_q1'; %Iq的坐标系在p上，再平移后得到Ip
Iq2 = R_p_to_q2 * principalmoments2 * R_p_to_q2'; %Iq的坐标系在p上，再平移后得到Ip
%
Ip1= Iq1 + m1*(p1'*p1*eye(3) - p1*p1');
Ip2= Iq2 + m2*(p2'*p2*eye(3) - p2*p2');
% [m1, m1*cx,m1*cy,m1*cz,Ixx1,Ixy1,Ixz1,Iyy1,Iyz1,Izz1]'
[m1;
m1*CoM1(1);
m1*CoM1(2);
m1*CoM1(3);
Ip1(1,1);
Ip1(1,2);
Ip1(1,3);
Ip1(2,2);
Ip1(2,3);
Ip1(3,3)
m2;
m2*CoM2(1);
m1*CoM2(2);
m1*CoM2(3);
Ip2(1,1);
Ip2(1,2);
Ip2(1,3);
Ip2(2,2);
Ip2(2,3);
Ip2(3,3)
]
%%
% Phi_1_cal = inv(A_1_1_num)*tau1_actual
% Phi_1_cal = inv(A_1_1_num'*A_1_1_num)*A_1_1_num'*tau1_actual
%自动svd分解方法，可行，
[U,S,V] = svd(A_1_1_num);
[r ,c] = size(S);
i=min(r,c);
i=3;
ApulsInverse = V*padarray(inv(S(1:i,1:i)),[c-i r-i],0,'post')*U';
Phi_1_cal = ApulsInverse*tau1_actual
%%
%手动svd分解方法，可行， m1*cx,m1*cz,Izz1
svdA_1_1_num = [A_1_1_num(:,2),A_1_1_num(:,3),A_1_1_num(:,10)];
Phi_1_cal = inv(svdA_1_1_num'*svdA_1_1_num)*svdA_1_1_num'*tau1_actual
Phi_1_cal = inv(A_1_1_num'*A_1_1_num)*A_1_1_num'*tau1_actual
%ridge regression岭回归方法，可行
Phi_1_cal = inv(A_1_1_num'*A_1_1_num + 0.0001*eye(10))*A_1_1_num'*tau1_actual% [m1, m1*cx,m1*cy,m1*cz,Ixx1,Ixy1,Ixz1,Iyy1,Iyz1,Izz1]'
tau1_cal = A_1_1_num * eval( [m1, m1*cx,m1*cy,m1*cz,Ixx1,Ixy1,Ixz1,Iyy1,Iyz1,Izz1]')

