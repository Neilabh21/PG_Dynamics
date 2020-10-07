% Prasanna Gandhi
% This is a code written to simulate cylinder applide by a force at one end
% cylinder parameter are given in the function file. 
% This program gives you animated look for
% the pendulum. 
% 
clf;
clear all;
l= 0.15;

% resolution 0 to 100s in 5000 steps
time = linspace(0,2,1000)';
% Initial Conditions  
IC = [0;0;0.010;0;0.010;0];
% Solve using ode45
[t,y] = ode45('Cylinderf',time,IC);
th = y(:,1); 
thd = y(:,2);
Xg = y(:,3);
Xgd = y(:,4);
Yg = y(:,5);
Ygd = y(:,6);

%%Plot of generalised coordinates w.r.t. time
figure(1)
plot(t,y(:,1)*180/(2*pi))
hold on
grid on
%plot(t,y(:,3),'b')
title('Cylinder angle');
xlabel('Time t');
ylabel('Angle (deg)');
legend('theta','r')
figure(2)
plot(y(:,3),y(:,5))
hold on
grid on
title('Cylinder trajectory')

%

%return
%% animation
figure(3)
O =[0 0];% origin or pivot point
axis(gca,'equal');% aspect ratio of the plot
axis([-0.1 0.3 -0.1 0.3]); % XY bounds
grid on;

% Loop for animation
for i = 1:length(t)
    % CG point X Y
    X= Xg(i); 
    Y= Yg(i);
    PCG = [X, Y];
    P1 = [X+l/2*cos(th(i)), Y+l/2*sin(th(i))];
    P2 = [X-l/2*cos(th(i)), Y-l/2*sin(th(i))];
    % circle at origin or pivot point
    %origincircle = viscircles(O,0.001);
    % Pendulum String joining pivot and current solution
    cylinder = line([P1(1) P2(1)],[P1(2) P2(2)], 'LineWidth',10);
     cg = line([PCG(1) PCG(1)],[PCG(2) PCG(2)], 'Marker','o','MarkerSize',2,'MarkerFaceColor','r');
    % ball
    %cg = viscircles(PCG, 0.007);
    % time interval to update the plot
    pause(0.05);
    %clear screen of previous pendulum position
    if i<length(t)
        %delete(origincircle);
        delete(cylinder);
        %delete(cg);
    end
end

