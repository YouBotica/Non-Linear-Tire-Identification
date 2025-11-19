% Save this as 'predict_NL_accel.m'
% This is the new "Grey Box" model
function [psi_ddot_nl] = predict_NL_accel(params, states, inputs, Fz_data, consts)
    
    % --- 1. Unpack Parameters (theta) ---
    % We now have 11 parameters
    Bf = params(1); Cf = params(2); Df = params(3); Ef = params(4); % Pacejka Front
    Br = params(5); Cr = params(6); Dr = params(7); Er = params(8); % Pacejka Rear
    P1 = params(9); P2 = params(10); P3 = params(11); % Residual Polynomial
    
    % --- 2. Unpack Constants ---
    lf = consts.lf;
    lr = consts.lr;
    Izz = consts.Izz;
    
    % --- 3. Unpack States (from smoothed data) ---
    vx = states(:, 1);
    vy = states(:, 2);
    r  = states(:, 3);
    
    % --- 4. Unpack Inputs ---
    delta = inputs(:, 2); 
    
    % --- 5. Get Vertical Loads (Fz) ---
    Fz_f_dyn = Fz_data(:, 1); 
    Fz_r_dyn = Fz_data(:, 2); 
    
    % --- 6. Calculate Slip Angles ---
    vx_safe = max(vx, 1.0);
    alpha_f = atan((vy + lf * r) ./ vx_safe) - delta;
    alpha_r = atan((vy - lr * r) ./ vx_safe);
    
    % --- 7. Calculate "White Box" (Pacejka) Forces ---
    F_peak_f = Df * Fz_f_dyn;
    x_f = Bf * alpha_f;
    Fyf = F_peak_f .* sin(Cf * atan(x_f - Ef * (x_f - atan(x_f))));
    
    F_peak_r = Dr * Fz_r_dyn;
    x_r = Br * alpha_r;
    Fyr = F_peak_r .* sin(Cr * atan(x_r - Er * (x_r - atan(x_r))));
    
    % --- 8. Calculate "White Box" (Pacejka) Acceleration ---
    psi_ddot_pacejka = (lf .* Fyf - lr .* Fyr) ./ Izz;
    
    % --- 9. Calculate "Black Box" (Residual) Acceleration ---
    % This part "soaks up" all the unmodeled physics (aero, etc.)
    % We model it as a polynomial of the states.
    psi_ddot_residual = (P1 * vx.^2) + (P2 * vy) + (P3 * r);

    % --- 10. Final Prediction ---
    % The total prediction is the sum of both models.
    psi_ddot_nl = psi_ddot_pacejka + psi_ddot_residual;
end