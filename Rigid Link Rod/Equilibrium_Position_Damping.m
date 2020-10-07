function [x, fval] = Equilibrium_Position_Damping()
% Equilibrium position of a beam
%   The rigid-rod torsion spring model is used to model the beam and
%   calculate the equilibrium position.

% Constants
t       = 0.4e-3;   % Thickness, m
w       = 0.02;     % Width, m
L       = 0.3;      % Length, m
g       = 9.81;     % Acceleration due to gravity, m/s^2
rho     = 8400;     % Density, kg/m^3
E       = 100e9;    % Young's Modulus

% Parameters
n_links = 3;      % No. of Links, dimensionless
P_by_P_cr = 10;     % Ratio of Mass to Critical Load

% Calculated Constants
m_rod   = rho * t * w * L;  % Mass of the rod, kg
l       = L / n_links;      % Length of each lnk, m
m       = m_rod / n_links;  % Mass of each link, kg
I       = t * w ^ 3 / 3;    % Moment of Inertia of cross section, m^4
k       = 3 * E * I / l;    % Spring constant for each link, N m / rad
P_cr    = pi ^2 * E * I / (4 * L ^ 2);  % Critical Load
M       = P_by_P_cr * P_cr / g;         % Tip mass

% Problem Definition
problem.options = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg', 'Display', 'iter-detailed'); 
% , 'Algorithm', 'trust-region-dogleg', 'Display', 'iter-detailed', 'PlotFcn', @optimplotfirstorderopt
% , 'Algorithm', 'levenberg-marquardt', 'Display', 'iter-detailed', 'PlotFcn', @optimplotfirstorderopt
% 'OptimalityTolerance', 1e-3, 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 4e3
problem.objective = @root2d;
problem.x0 = asin(0.3 * linspace(0, 1, n_links) .^ 2); % + 0.1 * rand([1, n_links]); % sin(linspace(0, 0.6, n_links) .* pi / l); % .^ 1; % +
problem.solver = 'fsolve';
[x,fval,exitflag,~] = fsolve(problem);
% disp(abs(x));
disp(['The state is ', num2str(x)]);
disp(['The function value is ', num2str(fval)]);
disp(['The exit flag is ', num2str(exitflag)]);
assert(exitflag == 1);

figure();
X0 = zeros([n_links+1,1]);
Y0 = zeros([n_links+1,1]);

for j=1:n_links
    X0(j+1) = X0(j) + l * sin(problem.x0(j));
    Y0(j+1) = Y0(j) + l * cos(problem.x0(j));
end

plot(X0, Y0);
xlim([-0.1, 0.1]);

figure();
X = zeros([n_links+1,1]);
Y = zeros([n_links+1,1]);

for j=1:n_links
    X(j+1) = X(j) + l * sin(x(j));
    Y(j+1) = Y(j) + l * cos(x(j));
end

plot(X, Y);
xlim([-0.1, 0.1]);


function F = root2d(x)

F(1) = k * (2 * x(1) - x(2)) - g * l * sin(x(1)) * (M + m * (n_links - 0.5));

if (n_links > 2)
    for i_1 = 2:(n_links-1)
        F(i_1) = k * (-x(i_1-1) + 2 * x(i_1) - x(i_1+1)) - g * l * sin(x(i_1)) * (M + m * (n_links - i_1 + 0.5));
    end
end

F(n_links) = k * (x(n_links) - x(n_links - 1)) - g * l * sin(x(n_links)) * (M + m * 0.5);

end
end
