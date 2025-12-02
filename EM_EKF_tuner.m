% --- 1. Load Data from Simulink Log ---
% Assumes you ran the sim and have 'out' in the workspace
if ~exist('out', 'var')
    error('Please run the Simulink simulation first to generate "out".');
end

% Vehicle Constants (Must match Simulink)
consts.m = veh_param.m;
consts.Izz = veh_param.Izz;
consts.lf = veh_param.lf;
consts.lr = veh_param.lr;
consts.Cd = veh_param.Cd;
% Pacejka Linear Slope (B*C*D)
consts.K_pacejka = veh_param.Bf * veh_param.Cf * veh_param.Df; 
dt = control_param.T;


% Run one pass with initial Q and R:
fprintf('Running Baseline (Human-Tuned) pass...\n');

% Extract Raw Measurements (Y) and Inputs (U)
% y = [vx_meas, vy_meas, r_meas, FzF_meas, FzR_meas]
% (IMPORTANT: Y should be measurements)
y_log = [out.measurements.Data(:,1), out.measurements.Data(:,2), out.measurements.Data(:,3), out.measurements.Data(:,4), out.FzF.Data, out.FzR.Data];

% u = [gp, delta]
u_log = out.u.Data;

Q_diag = control_param.ekf_Q_diag;

Q_curr = diag(Q_diag);

R_curr = diag([control_param.ekf_vx_sensor_noise; control_param.ekf_vy_sensor_noise; control_param.ekf_dpsi_sensor_noise; control_param.ekf_ay_sensor_noise]);

P_initial = control_param.ekf_initial_cov;


% Run one pass with your initial Q and R
[x_human_filt, P_human_filt, x_human_pred, P_human_pred, Ad_log] = ...
    forward_KF_EM(y_log, u_log, Q_curr, R_curr, P_initial, consts, dt);

% Run smoother on baseline
[x_human_smooth, ~, ~] = RTS_smoother_EM(x_human_filt, P_human_filt, x_human_pred, P_human_pred, Ad_log);

% Ground Truth for Validation (Yaw Accel)
gt_ddpsi = out.gt_ddpsi.Data;
gt_ddy = out.gt_ay.Data;

%%
[T_steps, ~] = size(y_log);

% --- 2. Initial Guess for Parameters ---
% Q: Process Noise Covariance
% We fix the top-left (physics states) to be small.
% We learn the bottom-right (disturbance states).
jitter = 1e-6;


% --- 3. EM Loop ---
MAX_ITER = 30; % 15-20 is usually enough for convergence
log_likelihood_history = zeros(MAX_ITER, 1);

fprintf('Starting EM Algorithm (%d iterations)...\n', MAX_ITER);

for iter = 1:MAX_ITER
    fprintf('  Iteration %d...', iter);
    
    % --- E-STEP: Forward Filter + RTS Smoother ---
    
    % 1. Run Forward EKF (MATLAB implementation)
    % fprintf(' Pre KF ');
    [x_filt, P_filt, x_pred, P_pred, Ad_log, log_likelihood] = forward_KF_EM(y_log, u_log, Q_curr, R_curr, control_param.ekf_initial_cov, consts, dt);
    % fprintf(' Post KF ');
    fprintf('    Log-Likelihood: %.4f\n', log_likelihood);


    % 2. Run RTS Smoother (with Lag-One Covariance)
    % fprintf(' Pre RTS ');
    [x_smooth, P_smooth, P_cross_smooth] = RTS_smoother_EM(x_filt, P_filt, x_pred, P_pred, Ad_log);
    % fprintf(' Post RTS ');
    
    % --- M-STEP: Update Q and R ---
    
    % 1. Update R (Measurement Noise)
   
    R_new_sum = zeros(4, 4);
    % fprintf(' Pre Rnew ');
    for t = 1:T_steps
        xt = x_smooth(t, :)'; % 6x1
        Pt = P_smooth(:, :, t); % 6x6
        yt = y_log(t, 1:4)'; % 3x1 (only the dynamic states)
        
        Fz = [out.FzF.Data(t); out.FzR.Data(t)];

        [Ad, Bd, Yv, Yr] = get_linear_model(xt, Fz, consts, dt);
        
        C_aug = [ % C matrix changes over time, more specifically its last row.
            eye(3), zeros(3);
            0, Yv, (Yr + xt(1)),  0, 1/consts.m, 0;    % Computes ay (NEW)
        ]; 
 

        % Residual
        res = yt - C_aug * xt;
        
        % Expectation term
        R_new_sum = R_new_sum + (res * res' + C_aug * Pt * C_aug');
    end
    % fprintf(' Post Rnew ');

    R_calculated = R_new_sum / T_steps;
    
    % 2. Update Q (Process Noise)
    Q_new_sum = zeros(6, 6);

    % fprintf(' Pre Qnew ');
    for t = 2:T_steps
        xt = x_smooth(t, :)';
        xt_prev = x_smooth(t-1, :)';
        ut_prev = u_log(t-1, :)';
        
        Pt = P_smooth(:, :, t);
        Pt_prev = P_smooth(:, :, t-1);
        P_cross = P_cross_smooth(:, :, t); % P_{t, t-1 | T}
        
        % Re-calculate A and B for time t-1 (Linearization)
        % [Ad, Bd, Yv, Yr] = get_linear_model(x, Fz, c, dt)
        [Ad, Bd, Yv, Yr] = get_linear_model(xt_prev, y_log(t-1,5:6), consts, dt);
        
        % Prediction residual (smoothed state vs model prediction)
        pred_res = xt - (Ad * xt_prev + Bd * ut_prev);
        
        % The "Correction" term involving covariances
        cov_correction = Pt + (Ad * Pt_prev * Ad') - (P_cross * Ad') - (Ad * P_cross');
        
        Q_new_sum = Q_new_sum + (pred_res * pred_res' + cov_correction);
    end
    % fprintf(' Post Qnew ');
    
    Q_calculated = Q_new_sum / (T_steps - 1);
    
    % --- CONSTRAINTS & UPDATES ---
    
    % 1. Force Diagonal structure
    Q_new = diag(diag(Q_calculated));
    R_new = diag(diag(R_calculated));
    
    % 2. Lock Physics States (Top-Left Block)?? NOTE: Maybe not
    % We trust our physics model for vx, vy, r.
    % Q_new(1,1) = control_param.ekf_Q_diag(1);
    % Q_new(2,2) = control_param.ekf_Q_diag(2);
    % Q_new(3,3) = control_param.ekf_Q_diag(3);
    
    % 3. Safety Floor
    % to prevent singularity in the next iteration.
    SAFE_FLOOR_Q = 1e-9; 
    SAFE_FLOOR_R = 1e-9;
    
    Q_new = max(Q_new, eye(6) * SAFE_FLOOR_Q);
    R_new = max(R_new, eye(4) * SAFE_FLOOR_R);
    
    % --- Update Parameters ---
    Q_curr = Q_new;
    R_curr = R_new; 
    
    % --- Validation Metric (RMSE of Lateral Accel Reconstruction) ---
    ddy_recon = zeros(T_steps, 1);
    
    for t = 1:T_steps
        % Unpack Smoothed States
        vx = max(x_smooth(t,1), 1.0);
        vy = x_smooth(t,2);
        r  = x_smooth(t,3);
        xi_y = x_smooth(t,5);
        
        % Unpack Inputs/Params
        delta = u_log(t,2);
        FzF = y_log(t,5); FzR = y_log(t,6); % Ensure indices match y_log structure
        
        % Recalculate dynamic stiffness (Total Axle Stiffness)
        Caf = consts.K_pacejka * FzF;
        Car = consts.K_pacejka * FzR;
        
        % --- 1. Linear Model (Body Frame dvy/dt) ---
        % Note: NO '2*' multiplier here!
        Yv = -(Caf + Car) / (consts.m * vx);
        Yr = -vx - (Caf*consts.lf - Car*consts.lr) / (consts.m * vx);
        B_lat = Caf / consts.m; 
        
        ay_linear_body = Yv*vy + Yr*r; % + B_lat*control_param.T*delta;
        
        % --- 2. Disturbance (Body Frame) ---
        ay_dist_body = xi_y / consts.m;
        
        % --- 3. Total Body Acceleration ---
        dvy_dt = ay_linear_body + ay_dist_body;
        
        % --- 4. Convert to IMU Acceleration (Inertial) ---
        % The IMU measures dvy/dt + vx*r (Centripetal term)
        % Note: Since Yr contains '-vx', adding '+vx*r' cancels it out,
        % leaving just the forces. This is what we want.
        ddy_recon(t) = dvy_dt + (vx * r);
    end
    
    % Calc RMSE against Ground Truth ay
    % (Make sure y_log(:,4) is actually ay in your log!)
    gt_ay = y_log(:, 4); 
    rmse = sqrt(mean((gt_ay - ddy_recon).^2));
    fprintf(' RMSE: %.4f\n', rmse);
    
end

fprintf('\nEM Algorithm Converged.\n');
disp('Learned Process Noise Q (Diagonal):');
disp(diag(Q_curr));
disp('Learned Measurement Noise R (Diagonal):');
disp(diag(R_curr));

% --- Plot Final Validation ---
figure;
plot(gt_ddy, 'r', 'LineWidth', 1.5); hold on;
plot(ddy_recon, 'b--', 'LineWidth', 1.5);
legend('Ground Truth', 'EM Reconstructed');
title('Validation: Reconstructed Yaw Acceleration after EM');
grid on;

fprintf('Plotting comparison...\n');
plot_EM_comparison(y_log, x_human_smooth, x_smooth, dt);

% --- HELPER: Linear Model Function ---
function [Ad, Bd, Yv, Yr] = get_linear_model(x, Fz, c, dt)
    vx = max(x(1), 1.0);
    FzF = Fz(1); FzR = Fz(2);

    Caf = c.K_pacejka * FzF;
    Car = c.K_pacejka * FzR;
    
    Yv = -(Caf + Car) / (c.m * vx);
    Yr = -vx - (Caf*c.lf - Car*c.lr) / (c.m * vx);
    Nv = -(Caf*c.lf - Car*c.lr) / (c.Izz * vx);
    Nr = -(Caf*c.lf^2 + Car*c.lr^2) / (c.Izz * vx);

    % disp('Yr');
    % disp(Yr);
    
    A = [ -c.Cd/c.m, 0, 0, 1/c.m, 0, 0;
          0, Yv, Yr, 0, 1/c.m, 0;
          0, Nv, Nr, 0, 0, 1/c.Izz;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0 ];
      
    B = [ 1/c.m, 0;
          0, Caf/c.m;
          0, Caf*c.lf/c.Izz;
          0, 0; 0, 0; 0, 0 ];
      
    Ad = eye(6) + A*dt;
    Bd = B*dt;
end