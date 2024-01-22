%%% This file creates an Mex file using the source code from C++ %%%

%% Clear all

clc; 
clear;

%% mex

%mex CXXFLAGS="$CXXFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" main.cpp economy.cpp initialization.cpp auxiliary.cpp 
mex CXXFLAGS="$CXXFLAGS -fopenmp -O3" LDFLAGS="$LDFLAGS -fopenmp" main.cpp economy.cpp initialization.cpp auxiliary.cpp
%mex main.cpp economy.cpp initialization.cpp auxiliary.cpp 

%% Common parameters:

params.b_grid_size_lowr = 99;          % Number of points in the grid for the low recovery bond. 
params.b_grid_size_highr = 101;         % Number of points in the grid for the high recovery bond. 
params.b_grid_min_lowr = 0.0;           % Minimum value for the low recovery bond grid.
params.b_grid_min_highr = 0.0;          % Minimum value for the high recovery bond grid.
params.b_grid_max_lowr = 1.2;           % Maximum value of the low recovery bond grid.
params.b_grid_max_highr = 1.2;          % Maximum value of the high recovery bond grid.
params.y_grid_size = 12;                 % Number of points in the grid for the income. (Always use Odd)
params.y_default = 0.969;               % Maximum income under default.
params.beta = 0.953;                    % Discount factor.
params.gamma = 2;                       % Risk aversion.
params.r = 0.017;                       % Interest rate.
params.rho = 0.945;                     % Persistence of the income.
params.sigma = 0.025;                   % Standard deviation of the income.
params.theta = 0.282;                   % Probability of a re-entry.
params.max_iter = 1000;                 % Maximum number of iterations.
params.tol = 1e-5;                      % Tolerance for the convergence.
params.m = 3;                           % Number of standard deviations for the income grid.
params.alpha_lowr = 0.00;               % Low recovery on defaulted debt.
params.alpha_highr = 0.15;              % High recovery on defaulted debt.

%% Run code with both alphas:

tic;
calibrated_model_solution = main(params);
save('Solution', 'calibrated_model_solution')
save('Parameters', 'params')
toc;