function dYdt = twolinkfun(t,Y)
k = 100;                   % Nm/radian
M = 0.1;                     % kg  
mbeam=0.1;                 % kg
m=mbeam/2;                 % kg
g = 9.81;                  % m/s^2
L=1;
l=L/2;
dYdt = zeros(2 , 1);
dYdt(1,1)=Y(2);
dYdt(2,1)=-(k*Y(1)+(m*g*l/2+M*g*l)*sin(Y(1)))/((M+m/3)*l^2);
end