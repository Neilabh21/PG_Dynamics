function x = one_iter()
% ================== %
% Defining Constants %
% ================== %

% ---------------- %
% Global Constants %
% ---------------- %

total_time = 0.01;                                                          % Seconds - Total Time
Del_t = 0.01;                                                               % Seconds - Time Step
Del_s = 0.01;                                                               % Metres - Element Size

N_t = total_time / Del_t + 1;                                               % Number of time steps 
N_s = 3;                                                                    % Number of spatial nodes
N_y = 9;                                                                    % Number of state variables

% --------- %
% Constants % 
% --------- %

E   = 0.05e9;                                                               % Pascals - Young's Modulus
% E = 210e9;                                                                  % Pascals - Young's Modulus
l   = 1;                                                                    % Metres - Length
nu  = 0.25;                                                                 % Unitless - Poisson Ratio
rho = 1500;                                                                 % kg/m^3 - Density
g   = 9.81;                                                                 % m/s^2 - Acceeleration due to gravity
R   = 0.005;                                                                % Metres - Radius
I   = pi * R^4 / 4;                                                         % Metres^4 - Moment of Inertia of the Cross Section
G   = E / (2 * (1 + nu));                                                   % * - Shear Modulus
k   = 0.9;                                                                  % Dimensionless - Roark Factor
A   = pi * R^2;                                                             % Metres^2 - Area
A_1 = k * A * G;                                                            % * - Shear Stiffness
A_2 = A * E;                                                                % * - Axial Stiffness
B   = E * I;                                                                % * - Bending Stiffness

% ------------------- %
% Smoothing Constants %
% ------------------- %

alpha_t = 0.5;
beta_t  = 0.5;
gamma_t = 0.5;
alpha_s = 0.5;
beta_s  = 0.5;
gamma_s = 0.5;

% --------------- %
% Define Matrices %
% --------------- %

M   = diag([1/A_1, 1/A_2, 1/(E * I), 0, 0, 0, rho * I, rho * A, rho * A]);  % M Matrix
K   = [ [0, 0, 0, 0, 0, 0, 0, 1, 0], ...
        [0, 0, 0, 0, 0, 0, 0, 0, 1], ...
        [0, 0, 0, 0, 0, 0, 1, 0, 0], ...
        [0, 0, 0, 1, 0, 0, 0, 0, 0], ...
        [0, 0, 0, 0, 1, 0, 0, 0, 0], ...
        [0, 0, 0, 0, 0, 1, 0, 0, 0], ...
        [0, 0, 1, 0, 0, 0, 0, 0, 0], ...
        [1, 0, 0, 0, 0, 0, 0, 0, 0], ...
        [0, 1, 0, 0, 0, 0, 0, 0, 0] ];                                      % K Matrix
M_c = diag([0, 0, 0, 1, 1, 1, 0, 0, 0]);                                    % M_c Matrix

% -------------- %
% Initialisation %
% -------------- %

y           = zeros(N_y, N_t, N_s);                                         % y - matrix
dely_delt   = zeros(N_y, N_t, N_s);                                         % del(y)/del(t) matrix
dely_dels   = zeros(N_y, N_t, N_s);                                         % del(y)/del(s) matrix
L           = zeros(N_y, N_t, N_s);                                         % Lambda Matrix
L_c         = zeros(N_y, N_t);                                              % Lambda_C Matrix

% =================== %
% Boundary conditions %
% =================== %

% Define later
f_1         = zeros(N_t, N_s);
f_2         = zeros(N_t, N_s);
m           = zeros(N_t, N_s);

% ========= %
% Main Loop %
% ========= %

options = optimoptions('fsolve','Display','off');
x0 = randn(4 * N_s + 1, N_y);
x0 = reshape(x0, (4 * N_s + 1) * N_y, 1);
for t_i = 2:N_t % t_i
    
    x = fsolve(@propagation, x0, options);
    % Update Y, L, L_c
    % calculate dely_delt, dely_dels
    % x = reshape(x, (2 * N_s + 1) * N_y);


end

% ===== %
% Plots %
% ===== %

% To be done later

% ========= %
% Functions %
% ========= %

function F =  propagation(x) % x : N_y x (4 * N_s + 1)     % at t_i
    
    x = reshape(x, 4 * N_s + 1, N_y);
    
    % --------------------- %
    % To Solve for 1st node %
    % --------------------- %
    
    % G^{i+1}_1
    F(N_s) = ((1 - alpha_t) / (gamma_t * Del_t)) * M_c * (x(:, 1) - y(:, t_i - 1, 1)) ...
                + ((alpha_t + gamma_t - 1) / (gamma_t * Del_t)) * M_c * dely_delt(:, t_i - 1, 1) ...
                + (1 - beta_t) * x(:, N_s + 1) ...
                + beta_t * y(:, t_i - 1, N_s + 1);
    
    % \Lambda^{i+1}_c
    F(N_s + 1) = x(:, N_s + 1) - [0, 0, 0, ...
                                    -x(7, 1), ...
                                    -x(8, 1) * cos(x(4, 1)) + x(9, 1) * sin(x(4,1)), ...
                                    -x(8, 1) * sin(x(4, 1)) * x(9, 1) * cos(x(4,1)), ...
                                    0, 0, 0]' ; % - Lambda_c in terms of x
    
    % ------------------------------------ %
    % To Solve for 2nd node to N_s^th node %
    % ------------------------------------ %
    
    for iter_S = 2:N_S
        F(iter_S-1) = M * ((1 - alpha_t) * ((1 - alpha_s) * x(:, 2 * N_s + 1 + iter_S) + alpha_s * x(:, 2 * N_s + iter_S)) ...
                            + alpha_t * ((1 - alpha_s) * dely_delt(:, t_i-1, iter_S) + alpha_s * dely_delt(:, t_i-1, iter_S-1))) ...
                    + K * ((1 - beta_t) * ((1 - beta_s) * x(:, 3 * N_s + 1 + iter_S) + beta_s * x(:, 3 * N_s + iter_S)) ...
                            + beta_t * ((1 - beta_s) * dely_dels(:, t_i-1, iter_S) + beta_s * dely_dels(:, t_i-1, iter_S-1))) ...
                    + ((1 - beta_t) * ((1 - beta_s) * x(:, 1 * N_s + 1 + iter_S) + beta_s * x(:, 1 * N_s + iter_S)) ...
                            + beta_t * ((1 - beta_s) * L(:, t_i-1, iter_S) + beta_s * L(:, t_i-1, iter_S-1)));
    end
    
    % -------------------- %
    % To Solve for Labda's % - To be updated
    % -------------------- %
    
    for iter_S = 1:N_s
        F(N_s + iter_S) = x(:, N_s + iter_S) - [-M/(E * I) * x(9, iter_S) + x(7, iter_S), ...
                                                M/(E * I) * x(8, iter_S), ...
                                                0, ...
                                                -M/(E * I), ...
                                                sin(x(4, iter_S)), ...
                                                -cos(x(4, iter_S)), ...
                                                -M/(E * I) * x(2, iter_S) + f_1(t_i, iter_S) + rho * A * x(7, iter_S) * x(9, iter_S), ...
                                                M/(E * I) * x(1, iter_S) + f_2(t_i, iter_S) - rho * A * x(7, iter_S) * x(8, iter_S), ...
                                                -x(2, iter_S) + m(t_i, iter_S)]';
    end
    
    % ---------------------- %
    % To Solve for dely_delt % ok
    % ---------------------- %
    
    for iter_S = 1:N_s
        F(2 * N_s + iter_S) = x(:, 2 * N_s + iter_S) ...
            - ((x(:, iter_S) - y(:, t_i-1, iter_S)) / (gamma_t * Del_t) - ((1 - gamma_t)/(gamma_t) * dely_delt(:, t_i-1, iter_S)));
    end
    
    % ---------------------- %
    % To Solve for dely_dels % kn
    % ---------------------- %
    
    F(3 * N_s + 1) = x(:, 3 * N_s + 1);
    
    for iter_S = 2:N_s
        F(3 * N_s + iter_S) = x(:, 3 * N_s + iter_S) ...
            - ((x(:, iter_S) - x(:, iter_S-1)) / (gamma_s * Del_s) - ((1 - gamma_s)/(gamma_s) * dely_dels(:, t_i-1, iter_S)));
    end
    
    % Reshaping to pass the function
    x = reshape(x, (4 * N_s + 1) * N_y, 1);
end
end
