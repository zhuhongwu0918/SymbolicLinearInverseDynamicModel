  function [c,ceq] = mycon(x)
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
   qF=qF+q{i};
   dqF=dqF+dq{i};
   ddqF=ddqF+ddq{i};
end

c=[abs(qF)-0.8;abs(dqF)-1.5];
% c=0;
ceq=0;

