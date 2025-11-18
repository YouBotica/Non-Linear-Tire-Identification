% Save this as 'predict_NL_accel.m'
function [psi_ddot_nl] = predict_NL_accel(params, states, inputs, consts)
    % This "Physics Whiz" model now correctly matches the
    % simulator's longitudinal load transfer physics.

    % --- 1. Unpack Parameters (theta) ---
    Bf = params(1); Cf = params(2); Df = params(3); Ef = params(4);
    Br = params(5); Cr = params(6); Dr = params(7); Er = params(8);
    
    % --- 2. Unpack Constants ---
    lf = consts.lf;
    lr = consts.lr;
    Izz = consts.Izz;
    m = consts.m;
    g = consts.g;
    h = consts.h;
    
    % --- 3. Unpack States (from smoothed data) ---
    vx = states(:, 1);
    vy = states(:, 2);
    r  = states(:, 3);
    
    % --- 4. Unpack Inputs ---
    gp = inputs(:, 1); % This is the longitudinal force Fx
    delta = inputs(:, 2);
    
    % --- 5. Calculate Slip Angles (Same) ---
    vx_safe = max(vx, 1.0);
    alpha_f = atan((vy + lf * r) ./ vx_safe) - delta;
    alpha_r = atan((vy - lr * r) ./ vx_safe);

    % --- 6. Calculate Vertical Loads (Fz) ---
    % This section now matches your simulator's equations
    
    % Calculate longitudinal acceleration
    % (We assume Fx is the dominant term, matching the sim)
    ax = gp / m; 
    
    % Calculate the longitudinal transfer term
    long_transfer_term = (ax - vy .* r) * m * h;
    
    % Calculate dynamic loads (ignoring aero for now)
    % This matches your F_zf and F_zr equations
    Fz_f_dyn = (lr * m * g - long_transfer_term) / (lf + lr);
    Fz_r_dyn = (lf * m * g + long_transfer_term) / (lf + lr);
    
    Fz_f_dyn = max(0, Fz_f_dyn);
    Fz_r_dyn = max(0, Fz_r_dyn);

    % --- 7. Calculate Tire Forces (Using the Fz-dependent Formula) ---
    % Front Tire
    F_peak_f = Df * Fz_f_dyn;
    x_f = Bf * alpha_f;
    Fyf = F_peak_f .* sin(Cf * atan(x_f - Ef * (x_f - atan(x_f))));
    
    % Rear Tire
    F_peak_r = Dr * Fz_r_dyn;
    x_r = Br * alpha_r;
    Fyr = F_peak_r .* sin(Cr * atan(x_r - Er * (x_r - atan(x_r))));
    
    % --- 8. Calculate Final Non-Linear Yaw Acceleration ---
    psi_ddot_nl = (lf .* Fyf - lr .* Fyr) ./ Izz;
end