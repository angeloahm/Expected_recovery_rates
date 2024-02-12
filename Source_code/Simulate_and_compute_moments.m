%% Compute the objects needed for the simulation:

Eq.B_policy_lowr = permute(reshape(calibrated_model_solution.B_policy_lowr + 1, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Eq.B_policy_highr = permute(reshape(calibrated_model_solution.B_policy_highr + 1, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Eq.D_policy = permute(reshape(calibrated_model_solution.D_policy, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);
Eq.B_grid_lowr = calibrated_model_solution.B_grid_lowr;
Eq.B_grid_highr = calibrated_model_solution.B_grid_highr;
Eq.Y_grid = calibrated_model_solution.Y_grid;
Eq.P = reshape(calibrated_model_solution.P, params.y_grid_size, params.y_grid_size)';


%% Simulate moments:
tic;
rng(0);
params_Simulation.TBurn = 10000;
params_Simulation.T = params_Simulation.TBurn + 1000*params_Simulation.TBurn;
Random_vec.theta = rand(params_Simulation.T,1);
[stats, simulated] = Run_Simulations(params, params_Simulation, Eq, Random_vec);
toc;