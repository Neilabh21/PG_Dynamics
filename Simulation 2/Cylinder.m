% Neilabh Banzal
% This is a code written to simulate cylinder applide by a force at one end
% cylinder parameter are given in the function file. 
% This program gives you animated look for
% the pendulum. 
% 
clf;
% clear all;
l= 0.15;

% resolution 0 to 100s in 5000 steps
time = linspace(0,2.5,1000)';
% Initial Conditions  
% IC = [0;0;0.010;0;0.010;0];
IC = zeros(3,1);
% Solve using ode45
[t,x_state] = ode45('Cylinderf',time,IC);
dt = (time(end) - time(1)) / length(time);
y_dot = x_state(:,1);
th = x_state(:,2); 
th_dot = x_state(:,3);

X_dot = - y_dot .* sin(th);
Y_dot = + y_dot .* cos(th);

X = zeros(length(time), 1);
Y = zeros(length(time), 1);

% X(1) = 10;
% Y(1) = -10;

for i = 2:length(time)
X(i) = X(i-1) - y_dot(i-1) * sin(th(i-1));
Y(i) = Y(i-1) + y_dot(i-1) * cos(th(i-1));
end

%%Plot of generalised coordinates w.r.t. time
% figure(1)
% plot(t, th * 180 / (pi))
% grid on
% %plot(t,y(:,3),'b')
% title('Cylinder angle');
% xlabel('Time t');
% ylabel('Angle (deg)');
% legend('theta','r')

figure(1)
% plot(t, X);
plot(X, Y);
grid on
title('Cylinder trajectory')
hold on
% plot(t, Y);
% legend('X', 'Y')

%

%return
%% animation
% figure(3)
% O =[0 0];% origin or pivot point
% axis(gca,'equal');% aspect ratio of the plot
% axis([-0.1 0.3 -0.1 0.3]); % XY bounds
% grid on;
% 
% % Loop for animation
% for i = 1:length(t)
%     % CG point X Y
%     X= Xg(i); 
%     Y= Yg(i);
%     PCG = [X, Y];
%     P1 = [X+l/2*cos(th(i)), Y+l/2*sin(th(i))];
%     P2 = [X-l/2*cos(th(i)), Y-l/2*sin(th(i))];
%     % circle at origin or pivot point
%     %origincircle = viscircles(O,0.001);
%     % Pendulum String joining pivot and current solution
%     cylinder = line([P1(1) P2(1)],[P1(2) P2(2)], 'LineWidth',10);
%      cg = line([PCG(1) PCG(1)],[PCG(2) PCG(2)], 'Marker','o','MarkerSize',2,'MarkerFaceColor','r');
%     % ball
%     %cg = viscircles(PCG, 0.007);
%     % time interval to update the plot
%     pause(0.05);
%     %clear screen of previous pendulum position
%     if i<length(t)
%         %delete(origincircle);
%         delete(cylinder);
%         %delete(cg);
%     end
% end
% 
