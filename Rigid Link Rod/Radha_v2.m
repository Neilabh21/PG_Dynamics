% initial conditions
n = 3;
rng(4);
Y0 = [deg2rad(10*(rand(n-1,1))); zeros(n-1,1)]; % Y contains 2 to n alphas and alpha_dots
tRange = linspace(0, 20, 40000);

% solving the ODE
[tSol, YSol] = ode45(@myODEfun_v2, tRange, Y0);

alpha=zeros((size(YSol,1)),n);
alpha(:,1)=zeros(size(alpha,1),1);
for i=2:n
    alpha(:,i)=YSol(:,i-1);
end

% coordinates of the bob in the base frame
l=1/n;
xb = zeros(size(YSol,1),1);
yb = zeros(size(YSol,1),1) + l;

for j=1: size(YSol,1)
for i=1:(size(YSol,2)/2)
    xb(j) = xb(j) + l * sin(YSol(j,i));
    yb(j) = yb(j) + l * cos(YSol(j,i));
end
end

% ball position

% trajectory plot
plot(xb,yb);


% animation
x = zeros(size(YSol,1), n+1); % x position of the start end of ith link
y = zeros(size(YSol,1), n+1); % y position of the start end of ith link
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

figure(2)
O =[0 0];% origin
axis(gca,'equal');% aspect ratio of the plot
axis([-1 1 -0.1 1.1]); % XY bounds
grid on;

for j = 1:100:length(tSol)-10
    xm=xb(j);  % x position of the pendulum bob
    ym=yb(j);  % y position of the pendulum bob
    bob_p=[xm, ym];
    
    for i=1:n
        rod(i)=line([x(j,i) x(j,(i+1))], [y(j,i) y(j,(i+1))], 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'LineWidth', 3);
    end
    hold on;
    plot(xb(1:j), yb(1:j), '-r');
    % ball
    %cg = viscircles(bob_p, 0.007);
    % time interval to update the plot
    pause(0.1);
    if j < length(tSol)
        for i=1:n
            delete(rod(i));
        end
    end
end



