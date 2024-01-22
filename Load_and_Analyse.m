%% Load results from previous computations:

load('Solution.mat')
load('Parameters.mat')

%% Format variables:

% We will reshape everything into cubes where:
%   - The x-axis is low recovery debt holdings (For example, from 0 to 0.6).
%   - The y-axis is high recovery debt holdings (For example, from 0 to 1.2).
%   - The z-axis is ouput level.

% Reshape exogenous objects:
Solution.Y_grid = calibrated_model_solution.Y_grid;
Solution.Y_grid_default = calibrated_model_solution.Y_grid_default;
Solution.B_grid_lowr = calibrated_model_solution.B_grid_lowr;
Solution.B_grid_highr = calibrated_model_solution.B_grid_highr;
Solution.P = reshape(calibrated_model_solution.P, params.y_grid_size, params.y_grid_size)';

% Included NaNs:
calibrated_model_solution.B_policy_lowr(calibrated_model_solution.D_policy == 1) = NaN;
calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy == 1) = NaN;

% Reshape endogenous objects:
Solution.Q_lowr = permute(reshape(calibrated_model_solution.Q_lowr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.Q_highr = permute(reshape(calibrated_model_solution.Q_highr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.V = permute(reshape(calibrated_model_solution.V, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.V_r = permute(reshape(calibrated_model_solution.V_r, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.V_d = permute(reshape(calibrated_model_solution.V_d, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.D_policy = permute(reshape(calibrated_model_solution.D_policy, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);

% When reshaping policy for debt, not that MATLAB starts indexing at 1.
Solution.B_policy_lowr = permute(reshape(calibrated_model_solution.B_policy_lowr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]) + 1;
Solution.B_policy_highr = permute(reshape(calibrated_model_solution.B_policy_highr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]) + 1;

%% Perform checks to see if bounds bind or not:

check_low_recovery = max(max(max(Solution.B_policy_lowr)));
check_high_recovery = max(max(max(Solution.B_policy_highr)));
disp(['Low recovery check: ' num2str(Solution.B_grid_highr(check_low_recovery))]);
disp(['High recovery check: ' num2str(Solution.B_grid_highr(check_high_recovery))]);



















%% Simulate moments:

rng(0);
X = 10000;
params_Simulation.TBurn = X;
params_Simulation.T = params_Simulation.TBurn + 1000*X;
Random_vec.theta = rand(params_Simulation.T,1);

[stats, simulated] = Run_Simulations(params, params_Simulation, Solution, Random_vec);

%% Plots 
hold on
plot(Solution.B_grid_lowr, Solution.B_policy_lowr(params.b_grid_size_highr,:,12))
plot(Solution.B_grid_highr, Solution.B_policy_highr(:,params.b_grid_size_lowr,12))