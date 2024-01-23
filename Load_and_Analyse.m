%% Reset:
clc 
clear all

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
disp(['Low recovery check: ' num2str(Solution.B_grid_lowr(check_low_recovery))]);
disp(['High recovery check: ' num2str(Solution.B_grid_highr(check_high_recovery))]);

%% Parameters for plots:

y_index = round(params.y_grid_size * 0.8);

bh_lowlevel = round(params.b_grid_size_highr/8);
bh_midlevel = round(params.b_grid_size_highr/4);
bh_highlevel = round(params.b_grid_size_highr/2);

bl_lowlevel = round(params.b_grid_size_lowr/8);
bl_midlevel = round(params.b_grid_size_lowr/4);
bl_highlevel = round(params.b_grid_size_lowr/2);


%% Plot value functions:


%%% Plot value functions for LOW_RECOVERY debt given high recovery debt:

% Set font size
fontSize = 14;
% Set figure size
figureSize = [10, 8];
% Set line width
lineWidth = 2;
% Set legend font size
legendFontSize = 12;
% Create figure
figure('Units', 'inches', 'Position', [0, 0, figureSize], 'Color', 'w');
% Create legend entry with LaTeX formatting
legend_entry_low = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_lowlevel));
legend_entry_mid = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_midlevel));
legend_entry_high = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_highlevel));
% Small high-recovery debt level:
plot(Solution.B_grid_lowr, Solution.V(bh_lowlevel,:, y_index), 'LineWidth', lineWidth, 'LineStyle', '-');
hold on;
% Medium high-recovery debt level:
plot(Solution.B_grid_lowr, Solution.V(bh_midlevel,:, y_index), 'LineWidth', lineWidth, 'LineStyle', '--');
hold on;
% Big high-recovery debt level:
plot(Solution.B_grid_lowr, Solution.V(bh_highlevel,:, y_index), 'LineWidth', lineWidth, 'LineStyle', '-.');
hold off;
% Add title and labels with LaTeX formatting
title('Value functions given $b_h$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$b_l$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('$V(\cdot, b_h, y)$', 'Interpreter', 'latex', 'FontSize', fontSize);
% Display the legend with increased font size
legend(legend_entry_low, legend_entry_mid, legend_entry_high, 'FontSize', legendFontSize);
% Display the grid
grid on;
% Increase font size of ticks
set(gca, 'FontSize', fontSize);
% Create folder if it doesn't exist
saveas(gcf, fullfile('Figures', 'value_functions_lowrecovery_plot.png'));

%%% Plot value functions for HIGH_RECOVERY debt given low recovery debt:

% Set font size
fontSize = 14;
% Set figure size
figureSize = [10, 8];
% Set line width
lineWidth = 2;
% Set legend font size
legendFontSize = 12;
% Create figure
figure('Units', 'inches', 'Position', [0, 0, figureSize], 'Color', 'w');
% Create legend entry with LaTeX formatting
legend_entry_low = sprintf('b_l = %.2f', Solution.B_grid_lowr(bl_lowlevel));
legend_entry_mid = sprintf('b_l = %.2f', Solution.B_grid_lowr(bl_midlevel));
legend_entry_high = sprintf('b_l = %.2f', Solution.B_grid_lowr(bl_highlevel));
% Small high-recovery debt level:
plot(Solution.B_grid_highr, Solution.V(:, bl_lowlevel, y_index), 'LineWidth', lineWidth, 'LineStyle', '-');
hold on;
% Medium high-recovery debt level:
plot(Solution.B_grid_highr, Solution.V(:, bl_midlevel, y_index), 'LineWidth', lineWidth, 'LineStyle', '--');
hold on;
% Big high-recovery debt level:
plot(Solution.B_grid_highr, Solution.V(:, bl_highlevel, y_index), 'LineWidth', lineWidth, 'LineStyle', '-.');
hold off;
% Add title and labels with LaTeX formatting
title('Value functions given $b_l$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$b_h$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('$V(b_l, \cdot, y)$', 'Interpreter', 'latex', 'FontSize', fontSize);
% Display the legend with increased font size
legend(legend_entry_low, legend_entry_mid, legend_entry_high, 'FontSize', legendFontSize);
% Display the grid
grid on;
% Increase font size of ticks
set(gca, 'FontSize', fontSize);
% Create folder if it doesn't exist
saveas(gcf, fullfile('Figures', 'value_functions_highrecovery_plot.png'));

%% Plot prices functions:

%%% Plot prices for LOW_RECOVERY debt given high recovery debt:

% Set font size
fontSize = 14;
% Set figure size
figureSize = [10, 8];
% Set line width
lineWidth = 2;
% Set legend font size
legendFontSize = 12;
% Create figure
figure('Units', 'inches', 'Position', [0, 0, figureSize], 'Color', 'w');
% Create legend entry with LaTeX formatting
legend_entry_low = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_lowlevel));
legend_entry_mid = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_midlevel));
legend_entry_high = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_highlevel));
% Small high-recovery debt level:
plot(Solution.B_grid_lowr, Solution.Q_lowr(bh_lowlevel,:, y_index), 'LineWidth', lineWidth, 'LineStyle', '-');
hold on;
% Medium high-recovery debt level:
plot(Solution.B_grid_lowr, Solution.Q_lowr(bh_midlevel,:, y_index), 'LineWidth', lineWidth, 'LineStyle', '--');
hold on;
% Big high-recovery debt level:
plot(Solution.B_grid_lowr, Solution.Q_lowr(bh_highlevel,:, y_index), 'LineWidth', lineWidth, 'LineStyle', '-.');
hold off;
% Add title and labels with LaTeX formatting
title('Low-recovery prices given $b_h$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$b_l$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('$Q_{l}(\cdot, b_h, y)$', 'Interpreter', 'latex', 'FontSize', fontSize);
% Display the legend with increased font size
legend(legend_entry_low, legend_entry_mid, legend_entry_high, 'FontSize', legendFontSize);
% Display the grid
grid on;
% Increase font size of ticks
set(gca, 'FontSize', fontSize);
% Create folder if it doesn't exist
saveas(gcf, fullfile('Figures', 'prices_lowrecovery_plot.png'));

%%% Plot prices for HIGH_RECOVERY debt given high recovery debt:

% Set font size
fontSize = 14;
% Set figure size
figureSize = [10, 8];
% Set line width
lineWidth = 2;
% Set legend font size
legendFontSize = 12;
% Create figure
figure('Units', 'inches', 'Position', [0, 0, figureSize], 'Color', 'w');
% Create legend entry with LaTeX formatting
legend_entry_low = sprintf('b_l = %.2f', Solution.B_grid_lowr(bl_lowlevel));
legend_entry_mid = sprintf('b_l = %.2f', Solution.B_grid_lowr(bl_midlevel));
legend_entry_high = sprintf('b_l = %.2f', Solution.B_grid_lowr(bl_highlevel));
% Small high-recovery debt level:
plot(Solution.B_grid_highr, Solution.Q_highr(:, bl_lowlevel, y_index), 'LineWidth', lineWidth, 'LineStyle', '-');
hold on;
% Medium high-recovery debt level:
plot(Solution.B_grid_highr, Solution.Q_highr(:, bl_midlevel, y_index), 'LineWidth', lineWidth, 'LineStyle', '--');
hold on;
% Big high-recovery debt level:
plot(Solution.B_grid_highr, Solution.Q_highr(:, bl_highlevel, y_index), 'LineWidth', lineWidth, 'LineStyle', '-.');
hold off;
% Add title and labels with LaTeX formatting
title('High recovery prices given $b_l$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$b_h$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('$Q_{h}(b_l, \cdot, y)$', 'Interpreter', 'latex', 'FontSize', fontSize);
% Display the legend with increased font size
legend(legend_entry_low, legend_entry_mid, legend_entry_high, 'FontSize', legendFontSize);
% Display the grid
grid on;
% Increase font size of ticks
set(gca, 'FontSize', fontSize);
% Create folder if it doesn't exist
saveas(gcf, fullfile('Figures', 'prices_highrecovery_plot.png.png'));

%% Plot policy functions:

%%% Plot policies for LOW_RECOVERY debt given high recovery debt:

% Set font size
fontSize = 14;
% Set figure size
figureSize = [10, 8];
% Set line width
lineWidth = 2;
% Set legend font size
legendFontSize = 12;
% Create figure
figure('Units', 'inches', 'Position', [0, 0, figureSize], 'Color', 'w');
% Create legend entry with LaTeX formatting
legend_entry_low = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_lowlevel));
legend_entry_mid = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_midlevel));
legend_entry_high = sprintf('b_h = %.2f', Solution.B_grid_highr(bh_highlevel));

% Small high-recovery debt level:
nonNaN_indices = ~isnan(Solution.B_policy_lowr(bh_lowlevel,:, y_index));
policy = Solution.B_policy_lowr(bh_lowlevel,:, y_index);
policy = policy(nonNaN_indices);
plot(Solution.B_grid_lowr(nonNaN_indices), Solution.B_grid_lowr(policy), 'LineWidth', lineWidth, 'LineStyle', '-');
hold on;
% Medium high-recovery debt level:
nonNaN_indices = ~isnan(Solution.B_policy_lowr(bh_midlevel,:, y_index));
policy = Solution.B_policy_lowr(bh_lowlevel,:, y_index);
policy = policy(nonNaN_indices);
plot(Solution.B_grid_lowr(nonNaN_indices), Solution.B_grid_lowr(policy), 'LineWidth', lineWidth, 'LineStyle', '--');
hold on;
% Big high-recovery debt level:
nonNaN_indices = ~isnan(Solution.B_policy_lowr(bh_highlevel,:, y_index));
policy = Solution.B_policy_lowr(bh_lowlevel,:, y_index);
policy = policy(nonNaN_indices);
plot(Solution.B_grid_lowr(nonNaN_indices), Solution.B_grid_lowr(policy), 'LineWidth', lineWidth, 'LineStyle', '-.');
hold off;
% Add title and labels with LaTeX formatting
title('Low-recovery bond policy given $b_h$', 'Interpreter', 'latex', 'FontSize', fontSize);
xlabel('$b_l$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('$b_{l}^{\prime}(\cdot, b_h, y)$', 'Interpreter', 'latex', 'FontSize', fontSize);
% Display the legend with increased font size
legend(legend_entry_low, legend_entry_mid, legend_entry_high, 'FontSize', legendFontSize);
% Display the grid
grid on;
% Increase font size of ticks
set(gca, 'FontSize', fontSize);
% Create folder if it doesn't exist
saveas(gcf, fullfile('Figures', 'policies_lowrecovery_plot.png'));



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