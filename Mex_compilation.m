%%% This file creates an Mex file using the source code from C++ %%%

%% Clear all
clc; 
clear;

%% mex
mex CXXFLAGS="$CXXFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" main.cpp economy.cpp initialization.cpp auxiliary.cpp 

%% Common parameters:
params.b_grid_size_lowr = 100;           % Number of points in the grid for the bond price.
params.b_grid_size_highr = 200;
params.b_grid_min_lowr = -0.8;         % Minimum value of the bond price.
params.b_grid_min_highr = -1.45;
params.b_grid_max_lowr = 0.0;         % Maximum value of the bond price.
params.b_grid_max_highr = 0.0;
params.y_grid_size = 21;               % Number of points in the grid for the income.
params.y_default = 0.969;              % Maximum income under default.
params.beta = 0.953;                   % Discount factor.
params.gamma = 2;                      % Risk aversion.
params.r = 0.017;                      % Interest rate.
params.rho = 0.945;                    % Persistence of the income.
params.sigma = 0.025;                  % Standard deviation of the income.
params.theta = 0.282;                  % Probability of a re-entry.
params.max_iter = 1000;                 % Maximum number of iterations.
params.tol = 1e-7;                     % Tolerance for the convergence.
params.m = 3;                          % Number of standard deviations for the income grid.
params.alpha_lowr = 0;                   % Low recovery on defaulted debt.
params.alpha_highr = 0.3;               % High recovery on defaulted debt.

%% Run code with both alphas;
tic;
calibrated_model_solution = main(params);
save('Solution', 'calibrated_model_solution')
save('Parameters', 'params')
toc;

%% Perform checks:

min(min(calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy==0)))
%% Format variables:

% Exogenous:
Y_grid = calibrated_model_solution.Y_grid;
Y_grid_default = calibrated_model_solution.Y_grid_default;
B_grid = calibrated_model_solution.B_grid;
P = reshape(calibrated_model_solution.P, params.y_grid_size, params.y_grid_size)';

%Endogenous: this stores everthing in 3-d matrices. Each matrix is:
%   |1, 2, 3, ..., nb  
%   |nb,nb+1,....,2nb-1,
%   |2nb,....

Q_low = permute(reshape(calibrated_model_solution.Q_low, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
Q_high = permute(reshape(calibrated_model_solution.Q_high, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
V = permute(reshape(calibrated_model_solution.V, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
V_r = permute(reshape(calibrated_model_solution.V_r, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
V_d = permute(reshape(calibrated_model_solution.V_d, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
B_policy_low = permute(reshape(calibrated_model_solution.B_policy_low, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
B_policy_high = permute(reshape(calibrated_model_solution.B_policy_high, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
D_policy = permute(reshape(calibrated_model_solution.D_policy, params.b_grid_size, params.b_grid_size, params.y_grid_size), [2, 1, 3]);
B_policy_high(D_policy == 1) = nan;
B_policy_low(D_policy==1) = nan;
B_policy_low = B_policy_low + 1;
B_policy_high = B_policy_high + 1;

%% Small analysis
figure(1)
plot(B_grid, B_policy_high(:,50,4));
figure(2)
plot(B_grid, D_policy(:,50,4));
figure(3)
plot(B_grid, D_policy(:,50,3));
figure(4)
plot(B_grid, Q_high(:,50,4));


