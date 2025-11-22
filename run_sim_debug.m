function run_sim_debug()
    clc; clear; close all;

    % =========================================================================
    % 1. CONFIGURATION & PARAMETERS
    % =========================================================================
    T = 0.01; % Sample time (100 Hz)
    SIM_TIME = 10.0; % Seconds
    N = floor(SIM_TIME / T);
    
    % --- Vehicle Parameters (780 kg Car) ---
    p.m   = 780;
    p.Izz = 1000;
    p.lf  = 1.4;
    p.lr  = 1.6;
    p.Cd  = 0.3; % Aero Drag Coeff
    p.g   = 9.81;
    p.h   = 0.35; % CG Height
    p.rho = 1.225;
    p.Af  = 1.5;  % Frontal Area
    
    % Pacejka Coeffs (Ground Truth)
    p.B = 22; p.C = 1.8; p.D = 1.6; p.E = 0.8;
    % Linear Stiffness Factor (for Controller/EKF)
    p.K_pacejka = p.B * p.C * p.D; 
    
    % --- Controller/EKF Parameters ---
    c.T = T;
    c.ekf_Q = diag([1e-6, 1e-6, 1e-6, 1000, 1000, 1000]); % Pure ADRC
    c.ekf_R = diag([0.001, 0.001, 0.001, 0.01]); % [vx, vy, r, ay]
    c.ekf_init_cov = eye(6);
    
    c.K1_vx = 4;
    c.K1_dpsi = 25;
    c.deccel_threshold = 5;
    
    % =========================================================================
    % 2. INITIALIZATION
    % =========================================================================
    % True State: [X, Y, psi, vx, vy, r]
    x_true = [0; 0; 0; 15; 0; 0]; 
    
    % EKF State: [vx, vy, r, xi_x, xi_y, xi_psi]
    x_hat = [15; 0; 0; 0; 0; 0];
    P = c.ekf_init_cov;
    
    % Inputs
    u = [0; 0]; % [gp, delta]
    
    % History Logging
    h.time = [];
    h.x_true = []; h.x_hat = []; 
    h.u = []; 
    h.ay_gt = []; h.ay_recon = [];
    
    fprintf('Starting Time-Domain Simulation (%d steps)...\n', N);
    
    % =========================================================================
    % 3. SIMULATION LOOP
    % =========================================================================
    for k = 1:N
        time = k*T;
        
        % --- A. REFERENCE GENERATION (Chirp) ---
        [psi_ref, vx_ref] = generate_references(time);
        
        % --- B. PLANT DYNAMICS (The "World") ---
        % 1. Apply inputs from previous step to get NEW state
        [x_next, measurements, debug_plant] = vehicle_plant_update(x_true, u, p, T);
        x_true = x_next;
        
        % 2. Add Sensor Noise
        noise = randn(4,1) .* [0.05; 0.01; 0.01; 0.1]; % vx, vy, r, ay
        y_meas = measurements + noise;
        
        % --- C. ESTIMATOR (EKF) ---
        % Note: We pass the Fz from the plant (perfect knowledge assumption for Fz 
        % matches your Simulink setup)
        Fz_meas = [debug_plant.FzF; debug_plant.FzR];
        
        [x_hat, P, debug_ekf] = run_EKF(x_hat, u, y_meas, Fz_meas, P, p, c);
        
        % --- D. CONTROLLER (ADRC) ---
        % Note: We need dpsi_ref for the inner loop. 
        % Simple approximation: dpsi_ref = K_heading * (psi_ref - psi)
        % Or if psi_ref is the chirp, its derivative is the dpsi_ref.
        % Let's use the P-controller approach you have:
        
        % Heading Error (Normalized)
        psi_err = angdiff(psi_ref, x_hat(3)); % Wait, x_hat(3) is r, not psi!
        % We need to integrate r to get psi for the controller? 
        % Or assumes x_hat tracks r and we use a separate heading.
        % For this DEBUG script, let's just control YAW RATE directly with the chirp.
        dpsi_ref_target = psi_ref; % (Re-purposing the variable for simplicity)
        
        [gp, delta] = run_ADRC(x_hat, vx_ref, dpsi_ref_target, p, c);
        u = [gp; delta];
        
        % --- E. LOGGING ---
        h.time(end+1) = time;
        h.x_true(end+1,:) = x_true';
        h.x_hat(end+1,:)  = x_hat';
        h.ay_gt(end+1)    = measurements(4); % True ay
        h.ay_recon(end+1) = debug_ekf.ay_recon;
        h.u(end+1,:)      = u';
    end
    
    % =========================================================================
    % 4. PLOTTING & DIAGNOSTICS
    % =========================================================================
    plot_results(h);
end

% =========================================================================
% HELPER FUNCTIONS (The "Guts")
% =========================================================================

function [x_next, y, dbg] = vehicle_plant_update(x, u, p, dt)
    % Unpack
    vx = x(4); vy = x(5); r = x(6);
    psi = x(3);
    gp = u(1); delta = u(2);
    
    % 1. Aero Forces
    Fx_aero = -0.5 * p.rho * p.Af * p.Cd * vx^2;
    
    % 2. Load Transfer (Longitudinal)
    ax_approx = gp / p.m; % Simplified for load transfer
    long_transfer = (ax_approx - vy*r) * p.m * p.h;
    
    FzF = (p.lr * p.m * p.g - long_transfer) / (p.lf + p.lr);
    FzR = (p.lf * p.m * p.g + long_transfer) / (p.lf + p.lr);
    
    % 3. Slip Angles (ISO Convention: delta - theta_v)
    % theta_v = atan2(vy + lr, vx)
    vx_safe = max(vx, 1.0);
    alpha_f = delta - atan2(vy + p.lf*r, vx_safe);
    alpha_r = 0     - atan2(vy - p.lr*r, vx_safe);
    
    % 4. Tire Forces (Pacejka '89)
    % Fy = D * sin(C * atan(B*alpha - E(B*alpha - atan(B*alpha))))
    FyF = p.D * FzF * sin(p.C * atan(p.B*alpha_f - p.E*(p.B*alpha_f - atan(p.B*alpha_f))));
    FyR = p.D * FzR * sin(p.C * atan(p.B*alpha_r - p.E*(p.B*alpha_r - atan(p.B*alpha_r))));
    
    % 5. Dynamics (Body Frame)
    % dvx = (Fx_tot / m) + vy*r
    % dvy = (Fy_tot / m) - vx*r
    % dr  = Mz / Izz
    
    Fx_total = gp + Fx_aero; % Assuming gp is net tractive force
    % Fy_total: Front projects onto Y axis
    Fy_total = FyF * cos(delta) + FyR; 
    
    dvx = (Fx_total / p.m) + vy*r;
    dvy = (Fy_total / p.m) - vx*r;
    
    % Moment
    % Front force creates moment: a * (FyF*cos(d))
    Mz = p.lf * (FyF * cos(delta)) - p.lr * FyR;
    dr = Mz / p.Izz;
    
    % 6. Kinematics (Global Frame)
    dX = vx * cos(psi) - vy * sin(psi);
    dY = vx * sin(psi) + vy * cos(psi);
    dpsi = r;
    
    % 7. Integration (Euler)
    x_next = x + [dX; dY; dpsi; dvx; dvy; dr] * dt;
    
    % 8. Measurements (The Sensors)
    % IMU ay = dvy + vx*r = (Fy_total/m - vx*r) + vx*r = Fy_total / m
    ay_true = Fy_total / p.m;
    
    y = [vx; vy; r; ay_true];
    
    dbg.FzF = FzF; dbg.FzR = FzR;
end

function [x_new, P_new, dbg] = run_EKF(x, u, y, Fz, P, p, c)
    % Unpack
    vx = x(1); vy = x(2); r = x(3);
    xi_x = x(4); xi_y = x(5); xi_psi = x(6);
    delta = u(2);
    
    % 1. Linear Model Update (Jacobians)
    vx_safe = max(vx, 1.0);
    
    % Calculate Stiffness from Fz inputs
    Caf = p.K_pacejka * Fz(1);
    Car = p.K_pacejka * Fz(2);
    
    Yv = -(Caf + Car) / (p.m * vx_safe);
    Yr = -vx - (Caf*p.lf - Car*p.lr) / (p.m * vx_safe); % Note the -vx
    Nv = -(Caf*p.lf - Car*p.lr) / (p.Izz * vx_safe);
    Nr = -(Caf*p.lf^2 + Car*p.lr^2) / (p.Izz * vx_safe);
    
    % A Matrix (Continuous)
    Ac = zeros(6,6);
    Ac(1,1) = -p.Cd/p.m; Ac(1,4) = 1/p.m;
    Ac(2,2) = Yv; Ac(2,3) = Yr; Ac(2,5) = 1/p.m;
    Ac(3,2) = Nv; Ac(3,3) = Nr; Ac(3,6) = 1/p.Izz;
    
    % B Matrix (Continuous)
    Bc = zeros(6,2);
    Bc(1,1) = 1/p.m;
    Bc(2,2) = Caf / p.m;
    Bc(3,2) = (Caf * p.lf) / p.Izz;
    
    % Discretize
    Ad = eye(6) + Ac * c.T;
    Bd = Bc * c.T;
    
    % 2. Predict
    x_pred = Ad * x + Bd * u;
    P_pred = Ad * P * Ad' + c.ekf_Q;
    
    % 3. Measurement Model (C Matrix)
    % y = [vx; vy; r; ay]
    % ay = Yv*vy + (Yr + vx)*r + xi_y/m + B(2,2)*delta (Wait, B term?)
    % Let's check ay reconstruction carefully:
    % dvy_model = Yv*vy + Yr*r + B(2,2)*delta + xi_y/m
    % ay_model = dvy_model + vx*r
    %          = Yv*vy + (Yr+vx)*r + B(2,2)*delta + xi_y/m
    
    C = zeros(4,6);
    C(1,1) = 1;
    C(2,2) = 1;
    C(3,3) = 1;
    C(4,2) = Yv; C(4,3) = (Yr + vx); C(4,5) = 1/p.m;
    
    % 4. Update
    % We must subtract the known input part B(2,2)*delta from the ay measurement
    % before innovation, OR include it in the prediction y_hat = C*x + D*u
    y_pred = C * x_pred;
    % Add the direct feedthrough term for ay (Steering -> Lateral Accel)
    ay_feedthrough = (Caf / p.m) * delta;
    y_pred(4) = y_pred(4) + ay_feedthrough;
    
    innovation = y - y_pred;
    
    S = C * P_pred * C' + c.ekf_R;
    K = P_pred * C' / S;
    
    x_new = x_pred + K * innovation;
    x_new(1) = max(x_new(1), 0.1); % Velocity floor
    P_new = (eye(6) - K*C) * P_pred;
    
    % 5. Debug Reconstruction
    % Reconstruct ay using the updated state
    % ay_recon = Yv*vy + (Yr+vx)*r + xi_y/m + ay_feedthrough
    dbg.ay_recon = C(4,:) * x_new + ay_feedthrough;
end

function [gp, delta] = run_ADRC(x, vx_ref, dpsi_ref, p, c)
    vx = x(1); vy = x(2); r = x(3);
    xi_x = x(4); xi_psi = x(6);
    
    % 1. Longitudinal
    U0_x = c.K1_vx * (vx_ref - vx) - (xi_x / p.m) + (p.Cd * vx / p.m);
    F_des = p.m * U0_x;
    if F_des > 0; gp = F_des;
    elseif F_des < -c.deccel_threshold; gp = F_des;
    else; gp = 0; end
    
    % 2. Lateral
    % Re-calc parameters
    vx_safe = max(vx, 1.0);
    % We assume nominal stiffness here for control to start
    % Ideally pass Fz in, but let's use static approximation for stability
    FzF_static = p.m*p.g*p.lr/(p.lf+p.lr);
    Caf = p.K_pacejka * FzF_static;
    Car = p.K_pacejka * (p.m*p.g - FzF_static);
    
    Nv = -(Caf*p.lf - Car*p.lr) / (p.Izz * vx_safe);
    Nr = -(Caf*p.lf^2 + Car*p.lr^2) / (p.Izz * vx_safe);
    B_lin = (Caf * p.lf) / p.Izz;
    
    U0_psi = c.K1_dpsi * (dpsi_ref - r) - (Nv*vy + Nr*r) - (xi_psi / p.Izz);
    delta = U0_psi / B_lin;
end

function [psi_ref, vx_ref] = generate_references(time)
    % Simple Sine Wave Yaw Rate
    vx_ref = 15;
    psi_ref = 0.5 * sin(2*pi*0.5*time); % 0.5 rad/s yaw rate target
end

function d = angdiff(a, b)
    d = a - b;
    d = mod(d + pi, 2*pi) - pi;
end

function plot_results(h)
    figure('Name', 'Simulation Debug', 'Color', 'w');
    
    subplot(3,1,1);
    plot(h.time, h.x_true(:,4), 'k--', 'LineWidth', 2); hold on;
    plot(h.time, h.x_hat(:,1), 'r-');
    ylabel('Vx'); legend('True', 'Est'); grid on;
    
    subplot(3,1,2);
    plot(h.time, h.x_true(:,6), 'k--', 'LineWidth', 2); hold on;
    plot(h.time, h.x_hat(:,3), 'r-');
    ylabel('Yaw Rate'); grid on;
    
    subplot(3,1,3);
    plot(h.time, h.ay_gt, 'k', 'LineWidth', 2); hold on;
    plot(h.time, h.ay_recon, 'c--');
    ylabel('Ay (Lateral Accel)');
    legend('Ground Truth', 'EKF Reconstructed');
    title('SANITY CHECK'); grid on;
end