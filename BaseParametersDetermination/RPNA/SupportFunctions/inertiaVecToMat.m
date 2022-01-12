function [ I ] = inertiaVecToMat( a )
%% [ I ] = inertiaVecToMat( a )
% Converts a vector of 10 inertial parameters to a spatial inertia matrix
% a=[m, hx, hy, hz, Ixx, Iyy, Izz, Iyz, Ixz, Ixy]';

Ibar = [a(5)  a(10) a(9) ; 
        a(10) a(6)  a(8) ; 
        a(9)  a(8)  a(7) ];%第二惯量
h    = [a(2)  a(3)  a(4)]';%第一惯量
m    = a(1);
I = [Ibar, skew(h) ; skew(h)', m*eye(3)];
%相对于坐标轴的刚体转动惯量的空间惯性矩阵，参考机器人手册中文版
%第2章 动力学，2.2.11空间惯量，33页
end

