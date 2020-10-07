% initial conditions
n = 3;
Y0 = [deg2rad(10*(rand(n-1,1))); zeros(n-1,1)];
tRange = [0 10];

% solving the ODE
[tSol, YSol] = ode45(@myODEfun, tRange, Y0);

% % global coordinates of the pendulum bob
% b = [0; 0 ;1];
% for i = 1:length(tSol)
%     rotation_matrix = reshape(YSol(i, 1:9), 3,3);
%     position = eul2rotm([0 0 pi]) * rotation_matrix * b;
%     Xcg(i) = position(1);
%     Ycg(i) = position(2);
%     Zcg(i) = position(3);
% end
% 
% 
% % animation
% figure(2)
% O =[0 0 0];% origin or pivot point
% axis(gca,'equal');% aspect ratio of the plot
% axis([-1 1 -1 1 -1 1]); % XYZ bounds
% grid on;
% 
% % Loop for animation
% for i = 1:length(tSol)
%     % CG point X Y Z
%     X = Xcg(i); 
%     Y = Ycg(i);
%     Z = Zcg(i);
%     PCG = [X, Y, Z];
%     % circle at origin or pivot point
%     % origincircle = viscircles(O,0.001);
%     % Pendulum String joining pivot and current solution
%     rod =   line([0 PCG(1)], [0 PCG(2)], [0 PCG(3)], 'Marker', '.', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'LineWidth', 3);
%     hold on;
%     plot3(Xcg(1:i), Ycg(1:i), Zcg(1:i), 'r');
%     % ball
%     % cg = viscircles(PCG, 0.007);
%     % time interval to update the plot
%     pause(0.1);
%     if i<length(tSol)
%         delete(rod);
%     end
% end

% function 
function dYdt = myODEfun(t,Y)
g = 9.81;                   % m/s^2
n_size = size(Y, 1) / 2;    % No units
l = 1 / (n_size + 1);       % Metres
k = 1;                      % Nm/radian
m = 1;                      % kg

dYdt = zeros(2 * n_size, 1);
for i = 1:n_size
    dYdt(i) = Y(n_size + i);
end

A = zeros(n_size, n_size);
for i = 1:n_size
    for j = 1:n_size
        A(i, j) = cos(Y(i) - Y(j));
    end
end

val_1 = 0;
val_2 = 0;
for i = 1:n_size
    val_1 = val_1 + Y(n_size + i) ^ 2 * sin(Y(i));
    val_2 = val_2 + Y(n_size + i) ^ 2 * cos(Y(i));
end

b = zeros(n_size, 1);
b(1) = - k / (m * l ^ 2) * (2 * Y(1) - 0 - Y(2)) + m * g / l * sin(Y(1)) + cos(Y(1)) * (val_1) - sin(Y(1)) * (val_2);
for i = 2:n_size-1
    b(i) = - k / (m * l ^ 2) * (2 * Y(i) - Y(i-1) - Y(i+1)) + m * g / l * sin(Y(i)) + cos(Y(i)) * (val_1) - sin(Y(i)) * (val_2);
end
b(n_size) = - k / (m * l ^ 2) * (2 * Y(n_size) - 0 - Y(n_size-1)) + m * g / l * sin(Y(n_size)) + cos(Y(n_size)) * (val_1) - sin(Y(n_size)) * (val_2);

dYdt(n_size + 1:end) = A\b;
end
