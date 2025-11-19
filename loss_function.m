% Save this as 'loss_function.m'

function [loss] = loss_function(params, states_smooth, inputs_log, Fz_data, consts, y_target, y_target_variance, initial_guess)
    
    % 1. Get the "Physics Whiz" prediction
    y_predict = predict_NL_accel(params, states_smooth, inputs_log, Fz_data, consts);
    
    % 2. Calculate the error
    error = y_target - y_predict;
    
    % --- 3. Calculate the *Weighted* Sum of Squared Errors ---
    safe_variance = max(y_target_variance, 1e-9);
    weighted_error_sq = (error.^2); %/ safe_variance;
    data_loss = sum(weighted_error_sq);

    % --- 4. (NEW) Calculate Regularization Penalty ---
    % We want to "pull" the parameters towards our initial_guess
    % (This is L2 Regularization)
    
    % We care more about the Pacejka params (1-8) than the residuals (9-11)
    % So we create a weight vector:
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, ... % Pacejka
               0.1, 0.1, 0.1];                           % Residuals (less penalty)
    
    lambda = 0.1; % <-- This is our "Penalty Strength" tuning knob
    
    param_error = (params - initial_guess) .* weights;
    regularization_loss = lambda * sum(param_error.^2);

    % --- 5. Final Loss ---
    % The total loss is the sum of the data error and the penalty
    loss = data_loss + regularization_loss;
    
    % Safety check
    if isnan(loss) || isinf(loss)
        loss = 1e100;
    end
end