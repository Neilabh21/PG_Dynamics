% constants
n = 2;
k = 100;                   % Nm/radian
M = 0.1;                     % kg  
mbeam=0.1;                 % kg
m=mbeam/n;                 % kg
g = 9.81;                  % m/s^2
l = 1 / n;                % Metres

% initial conditions
%rng(1);
%Y0 = [deg2rad(10*(rand(n-1,1)));deg2rad(10*(rand(n-1,1)));
%Y0 = [deg2rad(10*(rand(n-1,1))); ones(n-1,1)]; % Y contains 2 to n alphas and alpha_dots
%Y0 = [deg2rad(linspace(10,30,n-1)'); ones(n-1,1)]; % Y contains 2 to n alphas and alpha_dots
Y0=[deg2rad(30),10];
% Solving the ODE
%tRange =linspace(0,20,400000);
tRange=[0,10];
%tRange=[0,2.5];

%[tSol, YSol] = ode45(@myODEfun,tRange,Y0,M,m,n,g,k,l);
[tSol, YSol] = ode45(@myODEfunNewCorrected,tRange,Y0);

%% Alphas and Omegas
alpha=zeros((size(YSol,1)),n);
alpha(:,1)=zeros(size(alpha,1),1);
for i=2:n
    alpha(:,i)=YSol(:,i-1);
end

omegas=zeros((size(YSol,1)),n);
omegas(:,1)=zeros(size(omegas,1),1);
for i=2:n
    omegas(:,i)=YSol(:,n+i-2);
end

% coordinates of the bob in the base frame
l=1/n;
xb=zeros(size(YSol,1),1);
yb=l*ones(size(YSol,1),1);

for j=1: size(YSol,1)
for i=1:(size(YSol,2)/2)
    xb(j)=xb(j)+l*sin(YSol(j,i));
    yb(j)=yb(j)+l*cos(YSol(j,i));
end
end


% bob trajectory plot
figure(1);
plot(xb,yb);
title('bob trajectory');
xlabel('xb');
ylabel('yb');
hold off;



% alphas plot
figure(2);
for i=1:n-1
    plot(tSol,YSol(:,i));
    hold on;
end
title('Angles of all links with y axis');
xlabel('t');
ylabel('Alpha(rad)');
hold off;

%% energy plot
%Potential Energy
PE=zeros(length(tSol),1);
Pb=zeros(length(tSol),1); %bob
Ps=zeros(length(tSol),1); %spring
Pl=zeros(length(tSol),1); %link

%Kinetic energy
KE=zeros(length(tSol),1);
Kbc=zeros(length(tSol),1); %bob_cos term v
Kbs=zeros(length(tSol),1); %bob_sin term v
Klc=zeros(length(tSol),1); %link_cos term v
Kls=zeros(length(tSol),1); %link_sin term v
sq_omegas=zeros(length(tSol),1); %link omegas term

%Total Energy
TE=zeros(length(tSol),1);

for t=1:length(tSol)
    for i=2:n
        Pb(t)=Pb(t)+cos(alpha(t,i));
        if i==1
            Ps(t)=0;
        else
            Ps(t)=Ps(t)+(alpha(t,i)-alpha(t,i-1)).^2;
        end
        
       
        for j=1:i-1
            Pl(t)=Pl(t)+cos(alpha(t,j));
            Klc(t)=Klc(t)+omegas(t,j)*cos(alpha(t,j));
            Kls(t)=Kls(t)+omegas(t,j)*sin(alpha(t,j));
        end
        
        Pl(t)=Pl(t)+0.5*cos(alpha(t,i));
        Klc(t)=Klc(t)+0.5*omegas(t,i)*cos(alpha(t,i));
        Kls(t)=Kls(t)+0.5*omegas(t,i)*sin(alpha(t,i));
        
        Kbc(t)=Kbc(t)+omegas(t,i)*cos(alpha(t,i));
        Kbs(t)=Kbs(t)+omegas(t,i)*sin(alpha(t,i));
        
        sq_omegas(t)=sq_omegas(t)+(omegas(t,i)^2);
        
        
    end
    PE(t)=Pb(t)*M*g*l+Pl(t)*m*g*l+Ps(t)*k/2 + m*g*l/2;
    
    KE(t)=((M*l^2)/2)*(Kbc(t)^2+Kbs(t)^2) + ((m*l^2)/24)*sq_omegas(t) + ((m*l^2)/2)*(Klc(t)^2+Kls(t)^2);
    
    %Total Energy
    TE(t)=KE(t)+PE(t);
end


figure(3);
plot(tSol,TE);
title('Total Energy');
xlabel('t');
ylabel('TE');
hold off;

TE_percent_error=(max(TE)-min(TE))*100/mean(TE);
%% animation
x=zeros(size(YSol,1),n+1); % x position of the start end of ith link
y=zeros(size(YSol,1),n+1); % y position of the start end of ith link
%rod=line(n);
x(:,1)=zeros(size(YSol,1),1);
y(:,1)=zeros(size(YSol,1),1);
x(:,2)=zeros(size(YSol,1),1);
y(:,2)=l*ones(size(YSol,1),1);

x(:,n+1)=xb;
y(:,n+1)=yb;

for i=3:n
    x(:,i)=x(:,i-1)+l*sin(alpha(:,i));
    y(:,i)=y(:,i-1)+l*cos(alpha(:,i));
end

figure(4)
title('Animation & bob trajectory');
O =[0 0 0];% origin
axis(gca,'equal');% aspect ratio of the plot
axis([-1 1 -0.2 1.2]); % XYZ bounds
grid on;

for j = 1:100:length(tSol)
   xm=xb(j);  % x position of the pendulum bob
   ym=yb(j);  % y position of the pendulum bob
   zm=0;
   bob_p=[xm, ym, zm];
  
  for i=1:n
       rod(i)=line([x(j,i) x(j,(i+1))], [y(j,i) y(j,(i+1))], [0 0], 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'LineWidth', 3);
   end
   hold on;
    plot(xm, ym, 'r .');
    %ball
    %cg = viscircles(bob_p, 0.007);
    %time interval to update the plot
  pause(0.1);
   if j<length(tSol)
       for i=1:n 
           delete(rod(i));
       end
   end
end

hold off;
