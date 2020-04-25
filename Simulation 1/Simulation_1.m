% Defining Constants
F = 10;                                 % N
m = 2;                                  % kg
l = 0.2;                                % m
I = m * l * l / 12;                     % kg m^2
A = F * l / (4 * I);                    % s^-2

% Defining the system
syms t;                                 % Create cymbolic variable t
xint = int(sin(A * t^2), t, 0, t);      % Evaluating the first integral in symbolic notation
yint = int(cos(A * t^2), t, 0, t);      % Evaluating the first integral in symbolic notation

% Method 1                              % Evaluating the sybolic notation at each point
n = 10;                                 % Number of Iterations
x = zeros(n,1);                         % Defining the array to store the values
y = zeros(n,1);                         % Defining the array to store the values
for T = 1:n                             % Loop to calculate the values of the expression
x(T,1) = F * int(xint, t, 0, T) / m;
y(T,1) = F * int(yint, t, 0, T) / m;
end

% Method 2                              % Directly plotting from the symbolc equation
% Plotting x vs t
figure(1);
answerx = F * int(xint, t, 0, t) / m;   % Evaluating the second integral
fplot(answerx, [-2, 2]);                % Plotting the required function
title('$$x_t = \frac{F}{m} \int\displaylimits_0^t \left( \int\displaylimits_0^t sin(A t^2) dt \right) dt$$','interpreter','latex')
xlabel('t');
ylabel('x_t');
print('1. x vs t.jpg','-djpeg');

% Plotting y vs t
figure(2);
answery = F * int(yint, t, 0, t) / m;   % Evaluating the second integral
fplot(answery, [-2, 2]);                % Plotting the required function
title('$$y_t = \frac{F}{m} \int\displaylimits_0^t \left( \int\displaylimits_0^t cos(A t^2) dt \right) dt$$','interpreter','latex')
xlabel('t');
ylabel('y_t');
print('2. y vs t.jpg','-djpeg');

% Plotting x vs y
figure(3);
fplot(answerx, answery);
title('Trajectory');
xlabel('x');
ylabel('y');
print('3. x vs y.jpg','-djpeg');