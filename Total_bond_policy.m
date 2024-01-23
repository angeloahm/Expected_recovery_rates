calibrated_model_solution.B_policy_highr = calibrated_model_solution.B_policy_highr + 1;
calibrated_model_solution.B_policy_lowr = calibrated_model_solution.B_policy_lowr + 1;
Bond_policy = calibrated_model_solution.B_grid_highr(calibrated_model_solution.B_policy_highr) + calibrated_model_solution.B_grid_lowr(calibrated_model_solution.B_policy_lowr);
Bond_policy(calibrated_model_solution.D_policy == 1) = NaN;
Bond_policy = permute(reshape(Bond_policy, params.b_grid_size_lowr, params.b_grid_size_highr, params.y_grid_size), [2, 1, 3]);