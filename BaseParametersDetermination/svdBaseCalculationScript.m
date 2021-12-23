%%
%ref:M. Gautier, "Numerical calculation of the base inertial parameters of robots,"
% Journal of Robotic Systems, vol. 8, pp. 485-506, 1991.
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
Z = [2,0,2,5,0,4;9,0,2,5,0,18;...
    6,0,9,4,0,12;10,0,6,10,0,20;...
    1,0,6,4,0,2;5,0,2,2,0,10;...
    2,0,9,8,0,4; 10,0,7,4,0,20;...
    1,0,4,3,0,2;8,0,6,5,0,16;...
    9,0,5,1,0,18;9,0,1,2,0,18;...
    1,0,3,10,0,2; 4,0,2,10,0,8;...
    3,0,2,6,0,6;9,0,3,1,0,18;...
    5,0,5,3,0,10;10,0,1,4,0,20;...
    2,0,10,9,0,4;3,0,10,1,0,6];
%删除列为0的项
% Z = [2,2,4,5;9,2,18,5;...
%     6,9,12,4;10,6,20,10;...
%     1,6,2,4;5,2,10,2;...
%     2,9,4,8;10,7,20,4;...
%     1,4,2,3;8,6,16,5;...
%     9,5,18,1;9,1,18,2;...
%     1,3,2,10;4,2,8,10;...
%     3,2,6,6;9,3,18,1;...
%     5,5,10,3;10,1,20,4;...
%     2,10,4,9;3,10,6,1];
%%
%构造置换矩阵P
[Q_raw,R1_raw,~] = qr(Z);
[~,R1_diag] = qr(Z);
R1_diag = diag(R1_diag);
R1_diag = round(R1_diag,5);%可以使用精度阈值 size(Z,2)*max(R1_diag)*ksi,ksi是machine precision.
% R1_diag = numpy.linalg.qr(Z, mode='economic').diagonal().round(round)

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
dbi;
% fprintf('ddi =\n') 
ddi;
dbn =length(dbi);
n_dynparms=6%列数
P = eye(n_dynparms);%构造置换矩阵P
Pb = P(:,dbi);
Pd = P(:,ddi);
P = [Pb,Pd];
%%
[U,S,V] = svd(Z);
% S = thresholdsetzero(S,1e-4);
% V = thresholdsetzero(V,1e-4);
c=6;%列数
b=3;%特征值个数
r = 20;%行数
V1= V(:,1:3);%大小c*b 6*3
V2= V(:,4:end);%大小c*b 6*3
%Z*V = U*S
% Z*[V1,V2] = U*S
% Z*V2=0
PTV2= P'*V2;
V21 = PTV2(1:b,1:c-b);
V22 = PTV2(c-b+1:end,1:c-b);
beta = -V21*inv(V22);
syms as bs cs ds es fs real
dynparms = [as bs cs ds es fs]';
PTX = P'*dynparms;
X1 = PTX(1:3);
X2 = PTX(4:end);
beta = thresholdsetzero(beta,1e-4);
XB1 = X1 + beta * X2
