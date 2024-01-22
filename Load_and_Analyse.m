%% Load results from previous computations:

load('Solution.mat')
load('Parameters.mat')

%% Format variables:

Solution.Y_grid = calibrated_model_solution.Y_grid;
Solution.Y_grid_default = calibrated_model_solution.Y_grid_default;
Solution.B_grid_lowr = calibrated_model_solution.B_grid_lowr;
Solution.B_grid_highr = calibrated_model_solution.B_grid_highr;
Solution.P = reshape(calibrated_model_solution.P, params.y_grid_size, params.y_grid_size)';

calibrated_model_solution.B_policy_lowr(calibrated_model_solution.D_policy == 1 ) = NaN;
calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy == 1 ) = NaN;

% Endogenous objects:
Solution.Q_lowr = permute(reshape(calibrated_model_solution.Q_lowr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.Q_highr = permute(reshape(calibrated_model_solution.Q_highr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.V = permute(reshape(calibrated_model_solution.V, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.V_r = permute(reshape(calibrated_model_solution.V_r, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.V_d = permute(reshape(calibrated_model_solution.V_d, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Solution.B_policy_lowr = permute(reshape(calibrated_model_solution.B_policy_lowr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]) + 1;
Solution.B_policy_highr = permute(reshape(calibrated_model_solution.B_policy_highr, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]) + 1;
Solution.D_policy = permute(reshape(calibrated_model_solution.D_policy, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);

%% Take aways policies in default:


%% Perform checks to see if bounds bind or not:

lowest_high_check = min(calibrated_model_solution.B_policy_highr(calibrated_model_solution.D_policy==0));
lowest_low_check = min(calibrated_model_solution.B_policy_lowr(calibrated_model_solution.D_policy==0));
disp(['High recovery, Highest debt: ' num2str(calibrated_model_solution.B_grid_highr(lowest_high_check))]);
disp(['High recovery, Highest debt index: ' num2str(lowest_high_check)]);
disp(['Low recovery, Highest debt: ' num2str(calibrated_model_solution.B_grid_lowr(lowest_low_check))]);
disp(['Low recovery, Highest debt index:  ' num2str(lowest_low_check)]);

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