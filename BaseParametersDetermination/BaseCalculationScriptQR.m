%%
%ref:M. Gautier, "Numerical calculation of the base inertial parameters of robots," Journal of Robotic Systems, vol. 8, pp. 485-506, 1991.
clc
clear
close all
%%
%由回归矩阵K的符号表达式进行随机采样生成数值矩阵。
%生成矩阵a，第4列为0项，第5列为第一列*2,或直接使用给定的a矩阵
Z =  ceil(rand(20,3)*10);
Z(:,4) = zeros(20,1);
Z(:,5) = Z(:,1)*2;
% Z 20*6
% Z = [2,0,2,5,0,4;9,0,2,5,0,18;...
%     6,0,9,4,0,12;10,0,6,10,0,20;...
%     1,0,6,4,0,2;5,0,2,2,0,10;...
%     2,0,9,8,0,4; 10,0,7,4,0,20;...
%     1,0,4,3,0,2;8,0,6,5,0,16;...
%     9,0,5,1,0,18;9,0,1,2,0,18;...
%     1,0,3,10,0,2; 4,0,2,10,0,8;...
%     3,0,2,6,0,6;9,0,3,1,0,18;...
%     5,0,5,3,0,10;10,0,1,4,0,20;...
%     2,0,10,9,0,4;3,0,10,1,0,6];
%删除列为0的项
Z = [2,2,4,5;9,2,18,5;...
    6,9,12,4;10,6,20,10;...
    1,6,2,4;5,2,10,2;...
    2,9,4,8;10,7,20,4;...
    1,4,2,3;8,6,16,5;...
    9,5,18,1;9,1,18,2;...
    1,3,2,10;4,2,8,10;...
    3,2,6,6;9,3,18,1;...
    5,5,10,3;10,1,20,4;...
    2,10,4,9;3,10,6,1];
%%
%symbotics
syms as bs cs ds real
dynparms = [as bs cs ds]';
n_dynparms = size(dynparms,1);
[Q_raw,R1_raw,P_raw] = qr(Z);
[~,R1_diag] = qr(Z);
R1_diag = diag(R1_diag);
R1_diag = round(R1_diag,5);%可以使用精度阈值 size(Z,2)*max(R1_diag)*ksi,ksi是machine precision.

dbi=[];
ddi=[];
for i = 1:length(R1_diag)
    if R1_diag(i)~=0
        dbi = [dbi,i];
    else
        ddi=[ddi,i];
    end
end
% fprintf('dbi =\n') 
dbi
% fprintf('ddi =\n') 
ddi
dbn =length(dbi);

% dbi:R对角不为零的个数
P = eye(4);%构造置换矩阵P
Pb = P(:,dbi);
Pd = P(:,ddi);
P = [Pb,Pd];

[~,Rbd1] = qr(Z*P);%原矩阵列置换后的QR分解


Rb1 = Rbd1(1:dbn,1:dbn);%提取基
Rb1=round(Rb1,5);


Rd1 = Rbd1(1:dbn,dbn+1:end);
Rd1 = round(Rd1,5);


Kd = inv(Rb1)*Rd1;

base_idxs = [1:n_dynparms] * Pb; %提取基向量索引号

baseparms = round((Pb' + Kd*Pd'),5)*dynparms %生成一组最小惯性参数集
% as + 2*cs
%        bs
%        ds
n_base = length(baseparms)
