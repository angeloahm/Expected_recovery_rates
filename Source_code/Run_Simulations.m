function [stats, simulations] = Run_Simulations(p, p_sim, Eq, Rand_Vec)
    
    T = p_sim.T;
    TBurn = p_sim.TBurn;

    % Empty container:
    yt = NaN*zeros(T,1);
    yt_default = NaN*zeros(T,1);
    bt_lowr = NaN*zeros(T,1);
    bt_highr = NaN*zeros(T,1);
    d_t = NaN*zeros(T+1,1);

    % Generate the income path:
    current_y = round(p.y_grid_size/2); % Median income.
    current_b_lowr = 1;                 % No high recovery debt.
    current_b_highr = 1;                % No high recovery debt.
    d_t(1) = 0;

    for step = 1:T

        yt(step) = current_y;                                              % Current ouput.
        next_y = randsample(p.y_grid_size, 1, true, Eq.P(current_y, :));   % Stores de index for next period's output.

        if d_t(step) == 0      % No default.
            bt_lowr(step) = current_b_lowr;
            bt_highr(step) = current_b_highr;        
            next_bt_lowr = Eq.B_policy_lowr(current_b_highr, current_b_lowr, current_y);
            next_bt_highr = Eq.B_policy_highr(current_b_highr, current_b_lowr, current_y);
            d_t(step+1) = Eq.D_policy(next_bt_highr, next_bt_lowr, next_y); 
            current_b_lowr = next_bt_lowr;
            current_b_highr = next_bt_highr; 
            current_y = next_y;
        else                    % Default.
            if Rand_Vec.theta(step) > p.theta  % Fail to reenter.
                bt_lowr(step) = current_b_lowr;
                bt_highr(step) = current_b_highr;
                current_b_lowr = p.b_grid_size_lowr;
                current_b_highr = p.b_grid_size_highr;
                current_y = next_y;
                d_t(step+1) = 1;
            else
                bt_lowr(step) = current_b_lowr;
                bt_highr(step) = current_b_highr;
                d_t(step+1) = 0;
                current_b_lowr = 1; % All debt gets erased.
                current_b_highr = 1; % All debt gets erased.
                current_y = next_y;
            end
        end
    end

    simulations.Default_policy = d_t(TBurn:end-1);
    simulations.B_low = Eq.B_grid_lowr(bt_lowr(TBurn:end));
    simulations.B_low_index = bt_lowr(TBurn:end);
    simulations.B_high = Eq.B_grid_highr(bt_highr(TBurn:end));
    simulations.B_high_index = bt_highr(TBurn:end);
    simulations.B_total = simulations.B_low + simulations.B_high;
    simulations.Y = Eq.Y_grid(yt(TBurn:end));
    simulations.Y_index = yt(TBurn:end);

    simulations.B_highr_share = NaN(size(simulations.B_high));
    nonZeroTotalIndices = simulations.B_total ~= 0;
    simulations.B_highr_share(nonZeroTotalIndices) = simulations.B_high(nonZeroTotalIndices) ./ simulations.B_total(nonZeroTotalIndices);
    simulations.B_highr_share(simulations.Default_policy == 1) = NaN;

    stats.B_highr_share_mean = nanmean(simulations.B_highr_share);
    stats.B_highr_share_sd = nanstd(simulations.B_highr_share);

    stats.Y = mean(simulations.Y(simulations.Default_policy == 0));
    stats.B_lowr = mean(simulations.B_low(simulations.Default_policy == 0));
    stats.B_lowr_std = std(simulations.B_low(simulations.Default_policy == 0));
    stats.B_highr = mean(simulations.B_high(simulations.Default_policy == 0));
    stats.B_highr_std = std(simulations.B_high(simulations.Default_policy == 0));
    stats.B_total = mean(simulations.B_total(simulations.Default_policy == 0));
    stats.Default_policy = sum(simulations.Default_policy)/(T-TBurn);

end