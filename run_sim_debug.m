function run_sim_debug()
    clc; clear; close all;

    % =========================================================================
    % 1. CONFIGURATION
    % =========================================================================
    T = 0.01;       
    SIM_TIME = 5.0; 
    N = floor(SIM_TIME / T);
    
    % --- Vehicle Parameters ---
    p.m   = 780;    
    p.Izz = 1000;   
    p.lf  = 1.4;    
    p.lr  = 1.6;    
    p.Cd  = 0.3;    
    p.g   = 9.81;
    p.h   = 0.35;   
    p.rho = 1.225;
    p.Af  = 1.5;    

    % Pacejka Coeffs
    p.B = 22; p.C = 1.8; p.D = 1.6; p.E = 0.8;
    p.K_pacejka = p.B * p.C * p.D; 
    
    % --- Controller Parameters ---
    c.T = T;
    c.ekf_Q_diag = [1e-6; 1e-6; 1e-6; 1000; 1000; 1000]; 
    c.ekf_R_diag = [0.001; 0.001; 0.001; 0.01]; 
    c.ekf_init_cov = eye(6);
    
    c.K1_vx = 4;
    c.K1_dpsi = 25;
    c.deccel_threshold = 5;
    
    % =========================================================================
    % 2. INITIALIZATION (FIXED)
    % =========================================================================
    % True Plant State: [X, Y, psi, vx, vy, r] (6x1)
    % FIX: Must be 6 elements to match the plant function
    x_true = [0; 0; 0; 5; 0; 0]; 
    
    % EKF State: [vx, vy, r, xi_x, xi_y, xi_psi]
    x_hat = [5; 0; 0; 0; 0; 0]; 
    P = c.ekf_init_cov;
    
    u = [0; 0]; 
    
    % Logging
    log_time = zeros(N,1);
    log_ay_true = zeros(N,1);
    log_ay_recon = zeros(N,1);
    log_vx = zeros(N,1);
    
    fprintf('Starting Digital Twin Simulation...\n');

    % =========================================================================
    % 3. MAIN LOOP
    % =========================================================================
    for k = 1:N
        time = k*T;
        
        % --- A. Inputs ---
        if time < 1.0; delta_cmd = 0; else; delta_cmd = 0.05; end
        gp = 100; 
        u = [gp; delta_cmd];
        
        % --- B. PLANT DYNAMICS ---
        [x_next, measurements, dbg_plant] = vehicle_plant_step(x_true, u, p, T);
        x_true = x_next;
        
        % Add noise
        noise = (randn(4,1) .* sqrt(c.ekf_R_diag)); 
        y_meas = measurements + noise;
        
        % --- C. ESTIMATOR ---
        Fz_meas = [dbg_plant.FzF; dbg_plant.FzR];
        ekf_measurements = [y_meas; Fz_meas];
        
        [x_hat, P, ay_recon] = run_EKF_step(x_hat, u, ekf_measurements, P, p, c);
        
        % --- LOGGING ---
        log_time(k) = time;
        log_ay_true(k) = measurements(4); 
        log_ay_recon(k) = ay_recon;       
        log_vx(k) = x_true(4);
    end
    
    % =========================================================================
    % 4. PLOTTING
    % =========================================================================
    figure('Name', 'Digital Twin Diagnostics', 'Color', 'w');
    subplot(2,1,1);
    plot(log_time, log_ay_true, 'k', 'LineWidth', 2); hold on;
    plot(log_time, log_ay_recon, 'r--', 'LineWidth', 1.5);
    ylabel('Lateral Accel (m/s^2)');
    legend('Ground Truth', 'EKF Reconstruction'); grid on;
    
    subplot(2,1,2);
    plot(log_time, log_vx, 'b'); 
    ylabel('Vx (m/s)'); grid on;
end

% =========================================================================
% SUB-FUNCTIONS
% =========================================================================

function [x_next, y, dbg] = vehicle_plant_step(x, u, p, dt)
    vx = x(4); vy = x(5); r = x(6);
    psi = x(3);
    gp = u(1); delta = u(2);
    
    % 1. Loads
    ax = gp / p.m; 
    long_transfer = (ax - vy*r) * p.m * p.h;
    
    FzF = (p.lr * p.m * p.g - long_transfer) / (p.lf + p.lr);
    FzR = (p.lf * p.m * p.g + long_transfer) / (p.lf + p.lr);
    
    % FIX: Clamp Fz to be positive (Wheel Lift Protection)
    FzF = max(0, FzF);
    FzR = max(0, FzR);
    
    % 2. Slip Angles
    vx_safe = max(vx, 1.0);
    alpha_f = delta - atan2(vy + p.lf*r, vx_safe);
    alpha_r = 0     - atan2(vy - p.lr*r, vx_safe);
    
    % 3. Forces (Pacejka with ISO Sign: Positive Slip -> Negative Force)
    FyF = -p.D * FzF * sin(p.C * atan(p.B*alpha_f - p.E*(p.B*alpha_f - atan(p.B*alpha_f))));
    FyR = -p.D * FzR * sin(p.C * atan(p.B*alpha_r - p.E*(p.B*alpha_r - atan(p.B*alpha_r))));
    
    % 4. Dynamics
    FyF_body = FyF * cos(delta); 
    dvx = vy*r + gp/p.m;
    dvy = -vx*r + (FyF_body + FyR)/p.m;
    dr  = (p.lf * FyF_body - p.lr * FyR) / p.Izz;
    
    dX = vx * cos(psi) - vy * sin(psi);
    dY = vx * sin(psi) + vy * cos(psi);
    dpsi = r;
    
    x_next = x + [dX; dY; dpsi; dvx; dvy; dr] * dt;
    
    % 5. Sensor Output (ay = Fy_total / m)
    ay_sensor = (FyF_body + FyR) / p.m;
    
    y = [vx; vy; r; ay_sensor];
    dbg.FzF = FzF; dbg.FzR = FzR;
end

function [x_new, P_new, ay_recon] = run_EKF_step(x, u, meas_full, P, p, c)
    % Unpack
    meas = meas_full(1:4); 
    Fz   = meas_full(5:6);
    
    % 1. SAFE LINEARIZATION
    vx = max(x(1), 1.0); % Clamp velocity
    
    % Clamp Fz (prevent negative stiffness)
    Fz(1) = max(100, Fz(1)); 
    Fz(2) = max(100, Fz(2));
    
    Caf = p.K_pacejka * Fz(1);
    Car = p.K_pacejka * Fz(2);
    
    Yv = -(Caf + Car) / (p.m * vx);
    Yr = -vx - (Caf*p.lf - Car*p.lr) / (p.m * vx);
    Nv = -(Caf*p.lf - Car*p.lr) / (p.Izz * vx);
    Nr = -(Caf*p.lf^2 + Car*p.lr^2) / (p.Izz * vx);
    
    A = [ -p.Cd/p.m, 0, 0, 1/p.m, 0, 0;
          0, Yv, Yr, 0, 1/p.m, 0;
          0, Nv, Nr, 0, 0, 1/p.Izz;
          zeros(3,6) ];
          
    B = [ 1/p.m, 0; 0, Caf/p.m; 0, Caf*p.lf/p.Izz; zeros(3,2) ];
          
    Ad = eye(6) + A*c.T;
    Bd = B*c.T;
    
    % 2. PREDICT
    x_pred = Ad * x + Bd * u;
    
    % SAFE PRED: Check for NaN in x_pred
    if any(isnan(x_pred))
        warning('x_pred became NaN. Resetting to previous x.');
        x_pred = x; 
    end
    
    P_pred = Ad * P * Ad' + diag(c.ekf_Q_diag);
    
    % SAFE PRED: Symmetrize P_pred immediately
    P_pred = (P_pred + P_pred') / 2;

    % 3. MEASUREMENT MODEL
    C = zeros(4,6);
    C(1,1)=1; C(2,2)=1; C(3,3)=1;
    C(4,2)=Yv; C(4,3)=(Yr + vx); C(4,5)=1/p.m;
    
    % 4. UPDATE
    y_pred_state = C * x_pred;
    ay_feedthrough = (Caf/p.m) * u(2);
    y_pred_state(4) = y_pred_state(4) + ay_feedthrough;
    
    innovation = meas - y_pred_state;
    
    % SAFE INNOVATION: Clamp innovation to prevent massive jumps
    innovation = max(min(innovation, 10), -10); 
    
    % SAFE S Matrix: Add strong regularization
    S = (C * P_pred * C') + diag(c.ekf_R_diag) + eye(4)*1e-6;
    
    % Check condition number
    if rcond(S) < 1e-12
        warning('S matrix singular. Skipping update.');
        K = zeros(6,4);
    else
        K = P_pred * C' / S;
    end
    
    x_new = x_pred + K * innovation;
    
    % SAFE VELOCITY: Clamp output velocity
    x_new(1) = max(x_new(1), 0.1);
    
    I = eye(6);
    P_new = (I - K*C) * P_pred * (I - K*C)' + K*diag(c.ekf_R_diag)*K';
    P_new = (P_new + P_new') / 2; % Symmetrize
    
    % 5. RECONSTRUCT
    vy_new = x_new(2); r_new = x_new(3); xi_y_new = x_new(5);
    ay_recon = Yv*vy_new + (Yr+vx)*r_new + xi_y_new/p.m + ay_feedthrough;
end