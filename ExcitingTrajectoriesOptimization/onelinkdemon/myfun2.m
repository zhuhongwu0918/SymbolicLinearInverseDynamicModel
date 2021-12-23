function f=myfun2(x)
% f=x(1)^2+4*x(2)^2-3*x(2)*x(3)+0.5*x(3);
N=3;%num of Fourier
%  x=ones(2*N+1,1);
Tf=10;
omega=2*pi/Tf;
dt=0.007;
t=0:dt:Tf;

l=1;
q{1}=x(1)/(omega*l)*sin(omega*l*t)-x(2)/(omega*l)*cos(omega*l*t)+x(3);
dq{1}=x(1)*cos(omega*l*t)+x(2)*sin(omega*l*t);
ddq{1}=-x(1)*omega*l*sin(omega*l*t)+x(2)*omega*l*cos(omega*l*t);
for i=2:N
    l=i;
q{i}=x(2+(i-2)*2+2)/(omega*l)*sin(omega*l*t)-x(2+(i-2)*2+3)/(omega*l)*cos(omega*l*t);
dq{i}=x(2+(i-2)*2+2)*cos(omega*l*t)+x(2+(i-2)*2+3)*sin(omega*l*t);
ddq{i}=-x(2+(i-2)*2+2)*omega*l*sin(omega*l*t)+x(2+(i-2)*2+3)*omega*l*cos(omega*l*t);
end

qF=q{1};
dqF=dq{1};
ddqF=ddq{1};
for i=2:N
   qF=qF+q{i};%傅里叶级数叠加，获得q
   dqF=dqF+dq{i};%傅里叶级数叠加，获得dq
   ddqF=ddqF+ddq{i};%傅里叶级数叠加，获得ddq
end


F=[ddqF(1),dqF(1),sign(dqF(1)),cos(qF(1)),sin(qF(1))];
for i=2:length(t)
    F=[F;ddqF(i),dqF(i),sign(dqF(i)),cos(qF(i)),sin(qF(i))];%构造Y矩阵
end

%  x
% det(F'*F)
% f=-log(det(F'*F))
f=cond(F)






