clear ; clc; close all;
% 机器人各连杆SDH参数
alpha1=-0.5*pi;%
a1=0;
d1=0;
theta1=0;%变量

alpha2=-0.5*pi;
a2=0;
d2=185;%85
theta2=0;%变量

alpha3=0.5*pi;
a3=0;
d3=0;%变量
theta3=pi;

%           theta      d        a        alpha 
L(1)=Link([theta1      d1        a1        alpha1]); L(1).offset=0;L(1).qlim=[-pi,pi];
L(2)=Link([theta2      d2        a2        alpha2]); L(2).offset=-0.5*pi;L(2).qlim=[-pi,pi];%L(2).offset=-0.5*pi;
% 移动关节需要特别指定关节类型--jointtype
L(3)=Link([theta3      d3        a3        alpha3]);  L(3).qlim=[500,1000]; L(3).jointtype='P';
% 把上述连杆“串起来”
Scara=SerialLink(L,'name','LF');%'name','KirinLeftFrontLeg'%
% 定义机器人基坐标和工具坐标的变换
Scara.base = transl(0 ,0 ,0);
Scara.base = troty(0.5*pi);
% Scara.base = trotz(0.5*pi);
% Scara.base = trotz(0.5*pi);
Scara.gravity = [0,0,-9.8];
Scara.tool = transl(0 ,0 ,0);
Scara.teach();   
%给关节数值
% q=[0 0.5*pi 500 ];
% Scara.plot(q);
% 
% 
% q=[-pi/12 pi/6 500];
% Scara.fkine(q)
% 
% q2=Scara.ikine(ans)  %逆运动学
% j0=Scara.jacob0(q)    %雅可比矩阵