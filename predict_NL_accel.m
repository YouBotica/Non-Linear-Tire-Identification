% Save this as 'predict_NL_accel.m'
function [y_ddot_nl] = predict_NL_accel(params, states, inputs, Fz_data, consts)
    % This is the "Physics Whiz" model.
    % It now uses the Fz-dependent Pacejka formula to match the sim.

    % --- 1. Unpack Parameters (theta) ---
    Bf = params(1); Cf = params(2); Df = params(3); Ef = params(4);
    Br = params(5); Cr = params(6); Dr = params(7); Er = params(8);
    
    % --- 2. Unpack Constants ---
    lf = consts.lf;
    lr = consts.lr;
    Izz = consts.Izz;
    m = consts.m;
    
    % --- 3. Unpack States (from smoothed data) ---
    % x_hat = [vx, vy, r, xi_x, xi_y, xi_psi]'
    vx = states(:, 1);
    vy = states(:, 2);
    r  = states(:, 3);

    
    % --- 4. Unpack Inputs ---s
    delta = inputs(:, 2); % Assuming inputs_log = [gp, delta]

    % --- 5. Get Vertical Loads (Fz) ---
    % Get Fz values directly from the logged 'Fz_data' matrix
    Fz_f_dyn = Fz_data(:, 1); % Column 1 is FzF
    Fz_r_dyn = Fz_data(:, 2); % Column 2 is FzR

    % --- 6. Calculate Slip Angles --- NB: 
    vx_safe = max(vx, 1.0); % Prevent division by zero
    alpha_f = atan((vy + lf * r) ./ vx_safe) - delta; 
    alpha_r = atan((vy - lr * r) ./ vx_safe) - 0;

    % --- 7. Calculate Tire Forces (Using the *Correct* Fz-dependent Formula) ---  
    % Front Tire
    F_peak_f = Df * Fz_f_dyn;
    x_f = Bf * alpha_f;
    Fyf = -F_peak_f .* sin(Cf * atan(x_f - Ef * (x_f - atan(x_f))));
    
    % Rear Tire
    F_peak_r = Dr * Fz_r_dyn;
    x_r = Br * alpha_r;
    Fyr = -F_peak_r .* sin(Cr * atan(x_r - Er * (x_r - atan(x_r))));
    
    % --- 8. Calculate Final Non-Linear Yaw Acceleration ---
    % 
    y_ddot_nl = -vx.*r + (Fyf + Fyr) ./ m;
    psi_ddot_nl = (lf .* Fyf - lr .* Fyr) ./ Izz;
end