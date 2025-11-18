% Save this as 'loss_function.m'
function [loss] = loss_function(params, states, inputs, consts, y_target)
    % 1. Get the "Physics Whiz" prediction
    y_predict = predict_NL_accel(params, states, inputs, consts);
    
    % 2. Calculate the error
    error = y_target - y_predict;
    
    % 3. Calculate the sum of squared errors (our loss)
    loss = sum(error.^2);
end