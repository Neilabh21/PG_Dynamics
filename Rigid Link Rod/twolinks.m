k = 100;                   % Nm/radian
M = 0.1;                     % kg  
mbeam=0.1;                 % kg
m=mbeam/2;                 % kg
g = 9.81;                  % m/s^2
L=1;
l=L/2;

tRange=[0,10];
% initial conditions
Y0=[deg2rad(30),10];

[tsol,Y]=ode45(@twolinkfun,tRange,Y0);


xm=zeros(length(tsol),1);
ym=zeros(length(tsol),1);
for i=1:length(tsol)
    xm(i)=l*sin(Y(i,1));
    ym(i)=l*cos(Y(i,1));
end
ym_=l*ones(length(tsol),1)+ym;
%energy
TE=zeros(length(tsol),1);
for t=1:length(tsol)
TE(t,1)=(m*g*l/2)+(m*g*l+m*g*l*cos(Y(t,1))/2)+(0.5/3)*m*l^2*Y(t,2)^2+(0.5*M*l^2*Y(t,2)^2)+M*g*ym_(t);
end

plot(tsol,TE);
title('Total Energy');
xlabel('t');
ylabel('TE');
hold off;

figure (2) 
plot(xm,ym_);
    