function [loss] = loss_function_weighted(params, states_smooth, inputs_log, Fz_data, consts, y_target, y_target_variance, initial_guess)
    
    % --- 1. Get the "Physics Whiz" Prediction ---
    y_predict = predict_NL_accel(params, states_smooth, inputs_log, Fz_data, consts);
    
    % --- 2. Calculate Weighted Data Loss (Likelihood) ---
    residuals = y_target - y_predict;
    safe_variance = max(y_target_variance, 1e-6); 
    data_loss = sum((residuals.^2) ./ safe_variance);
    
    % --- 3. Calculate Prior Regularization (Safety Net) ---
    % "Don't stray too far from the initial guess"
    lambda_prior = 0.1; % Keep this relatively loose
    
    scale_factors = abs(initial_guess);
    scale_factors(scale_factors < 1e-6) = 1.0;
    
    % Weights: Punish residuals (9-11) less than Pacejka (1-8)
    prior_weights = [1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1];
    
    param_dev = ((params - initial_guess) ./ scale_factors) .* prior_weights;
    prior_loss = lambda_prior * sum(param_dev.^2);

    % --- 4. (NEW) Calculate Similarity Loss (Front vs Rear) ---
    % "Front tires should behave similarly to Rear tires"
    
    lambda_sim = 100; % <-- TUNING KNOB: Higher = Force them closer together
    
    % Extract Pacejka params
    p_front = params(1:4); % [Bf, Cf, Df, Ef]
    p_rear  = params(5:8); % [Br, Cr, Dr, Er]
    
    % Use initial guess scaling so B (20) and D (1.0) are treated equally
    sim_scale = abs(initial_guess(1:4));
    
    % Calculate the difference
    sim_diff = (p_front - p_rear) ./ sim_scale;
    
    similarity_loss = lambda_sim * sum(sim_diff.^2);

    % --- 5. Total Loss ---
    loss = data_loss + prior_loss + similarity_loss;
    
    % Safety
    if isnan(loss) || isinf(loss), loss = 1e20; end
end