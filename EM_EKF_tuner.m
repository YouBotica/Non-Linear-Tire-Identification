% --- 1. Load Data from Simulink Log ---
% Assumes you ran the sim and have 'out' in the workspace
if ~exist('out', 'var')
    error('Please run the Simulink simulation first to generate "out".');
end

% Vehicle Constants (Must match Simulink)
consts.m = 780;
consts.Izz = 1000;
consts.lf = 1.4;
consts.lr = 1.6;
consts.Cd = 0.3;
% Pacejka Linear Slope (B*C*D)
consts.K_pacejka = 22 * 1.8 * 1.6; 
dt = 0.01;


% Run one pass with initial Q and R:
fprintf('Running Baseline (Human-Tuned) pass...\n');

% Extract Raw Measurements (Y) and Inputs (U)
% y = [vx_meas, vy_meas, r_meas, FzF_meas, FzR_meas]
% (IMPORTANT: Y should be measurements)
y_log = [out.measurements.Data(:,1), out.measurements.Data(:,2), out.measurements.Data(:,3), out.FzF.Data, out.FzR.Data];

% u = [gp, delta]
u_log = out.u.Data;

Q_diag = [
    1e-9;  % vx (Trust model)
    1e-6;  % vy (Trust model)
    1e-6;  % r (Trust model)
    10000;  % xi_x (Disturbance can change a lot!)
    10000;  % xi_y (Disturbance can change a lot!)
    1000; % xi_psi (Disturbance can change a lot!)
];

Q_curr = diag(Q_diag);

R_curr = diag([0.01; 0.01; 0.001]);


% Run one pass with your initial Q and R
[x_human_filt, P_human_filt, x_human_pred, P_human_pred, Ad_log] = ...
    forward_KF_EM(y_log, u_log, Q_curr, R_curr, consts, dt);

% Run smoother on baseline
[x_human_smooth, ~, ~] = RTS_smoother_EM(x_human_filt, P_human_filt, x_human_pred, P_human_pred, Ad_log);


% Ground Truth for Validation (Yaw Accel)
gt_ddpsi = out.gt_ddpsi.Data;

[T_steps, ~] = size(y_log);
dt = 0.01;

% --- 2. Initial Guess for Parameters ---
% Q: Process Noise Covariance
% We fix the top-left (physics states) to be small.
% We learn the bottom-right (disturbance states).
jitter = 1e-9;
Q_curr = diag([jitter, jitter, jitter, 10000, 100, 100]);

% R: Measurement Noise Covariance
% Initial guess: variance = 0.01 (std dev = 0.1)
R_curr = diag([0.01, 0.01, 0.01]);


% --- 3. EM Loop ---
MAX_ITER = 20; % 15-20 is usually enough for convergence
log_likelihood_history = zeros(MAX_ITER, 1);

fprintf('Starting EM Algorithm (%d iterations)...\n', MAX_ITER);

for iter = 1:MAX_ITER
    fprintf('  Iteration %d...', iter);
    
    % --- E-STEP: Forward Filter + RTS Smoother ---
    
    % 1. Run Forward EKF (MATLAB implementation)
    [x_filt, P_filt, x_pred, P_pred, Ad_log] = forward_KF_EM(y_log, u_log, Q_curr, R_curr, consts, dt);
    
    % 2. Run RTS Smoother (with Lag-One Covariance)
    [x_smooth, P_smooth, P_cross_smooth] = RTS_smoother_EM(x_filt, P_filt, x_pred, P_pred, Ad_log);
    
    % --- M-STEP: Update Q and R ---
    
    % 1. Update R (Measurement Noise)
    C = [eye(3), zeros(3,3)]; 
    R_new_sum = zeros(3, 3);
    
    for t = 1:T_steps
        xt = x_smooth(t, :)'; % 6x1
        Pt = P_smooth(:, :, t); % 6x6
        yt = y_log(t, 1:3)'; % 3x1 (only the dynamic states)
        
        % Residual
        res = yt - C * xt;
        
        % Expectation term
        R_new_sum = R_new_sum + (res * res' + C * Pt * C');
    end
    R_calculated = R_new_sum / T_steps;
    
    % 2. Update Q (Process Noise)
    Q_new_sum = zeros(6, 6);
    
    for t = 2:T_steps
        xt = x_smooth(t, :)';
        xt_prev = x_smooth(t-1, :)';
        ut_prev = u_log(t-1, :)';
        
        Pt = P_smooth(:, :, t);
        Pt_prev = P_smooth(:, :, t-1);
        P_cross = P_cross_smooth(:, :, t); % P_{t, t-1 | T}
        
        % Re-calculate A and B for time t-1 (Linearization)
        [Ad, Bd] = get_linear_model(xt_prev, ut_prev, u_log(t-1,2), y_log(t-1,4:5), consts, dt);
        
        % Prediction residual (smoothed state vs model prediction)
        pred_res = xt - (Ad * xt_prev + Bd * ut_prev);
        
        % The "Correction" term involving covariances
        cov_correction = Pt + (Ad * Pt_prev * Ad') - (P_cross * Ad') - (Ad * P_cross');
        
        Q_new_sum = Q_new_sum + (pred_res * pred_res' + cov_correction);
    end
    
    Q_calculated = Q_new_sum / (T_steps - 1);
    
    % --- CONSTRAINTS & UPDATES ---
    
    % 1. Force Diagonal structure
    Q_new = diag(diag(Q_calculated));
    R_new = diag(diag(R_calculated));
    
    % 2. Lock Physics States (Top-Left Block)
    % We trust our physics model for vx, vy, r.
    Q_new(1,1) = jitter;
    Q_new(2,2) = jitter;
    Q_new(3,3) = jitter;
    
    % 3. Safety Floor (THE CRITICAL FIX)
    % Clamp Q and R *before* assigning to Q_curr/R_curr
    % to prevent singularity in the next iteration.
    SAFE_FLOOR_Q = 1e-6; 
    SAFE_FLOOR_R = 1e-6;
    
    Q_new = max(Q_new, eye(6) * SAFE_FLOOR_Q);
    R_new = max(R_new, eye(3) * SAFE_FLOOR_R);
    
    % --- Update Parameters ---
    Q_curr = Q_new;
    R_curr = R_new;
    
    % --- Validation Metric (RMSE of Yaw Accel Reconstruction) ---
    ddpsi_recon = zeros(T_steps, 1);
    for t = 1:T_steps
        vx = max(x_smooth(t,1), 1.0);
        vy = x_smooth(t,2);
        r  = x_smooth(t,3);
        xi_psi = x_smooth(t,6);
        delta = u_log(t,2);
        FzF = y_log(t,4); FzR = y_log(t,5);
        
        % Recalculate dynamic stiffness
        Caf = consts.K_pacejka * FzF;
        Car = consts.K_pacejka * FzR;
        
        Nv = -(2*Caf*consts.lf - 2*Car*consts.lr) / (consts.Izz * vx);
        Nr = -(2*Caf*consts.lf^2 + 2*Car*consts.lr^2) / (consts.Izz * vx);
        B_lin = (2*Caf*consts.lf) / consts.Izz;
        
        ddpsi_recon(t) = (Nv*vy + Nr*r + B_lin*delta) + (xi_psi / consts.Izz);
    end
    
    rmse = sqrt(mean((gt_ddpsi - ddpsi_recon).^2));
    fprintf(' RMSE: %.4f\n', rmse);
    
end

fprintf('\nEM Algorithm Converged.\n');
disp('Learned Process Noise Q (Diagonal):');
disp(diag(Q_curr));
disp('Learned Measurement Noise R (Diagonal):');
disp(diag(R_curr));

% --- Plot Final Validation ---
figure;
plot(gt_ddpsi, 'r', 'LineWidth', 1.5); hold on;
plot(ddpsi_recon, 'b--', 'LineWidth', 1.5);
legend('Ground Truth', 'EM Reconstructed');
title('Validation: Reconstructed Yaw Acceleration after EM');
grid on;

fprintf('Plotting comparison...\n');
plot_EM_comparison(y_log, x_human_smooth, x_smooth, dt);

% --- HELPER: Linear Model Function ---
function [Ad, Bd] = get_linear_model(x, u, delta, Fz, c, dt)
    vx = max(x(1), 1.0);
    FzF = Fz(1); FzR = Fz(2);
    
    Caf = c.K_pacejka * FzF;
    Car = c.K_pacejka * FzR;
    
    Yv = -(2*Caf + 2*Car) / (c.m * vx);
    Yr = -vx - (2*Caf*c.lf - 2*Car*c.lr) / (c.m * vx);
    Nv = -(2*Caf*c.lf - 2*Car*c.lr) / (c.Izz * vx);
    Nr = -(2*Caf*c.lf^2 + 2*Car*c.lr^2) / (c.Izz * vx);
    
    A = [ -c.Cd/c.m, 0, 0, 1/c.m, 0, 0;
          0, Yv, Yr, 0, 1/c.m, 0;
          0, Nv, Nr, 0, 0, 1/c.Izz;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0 ];
      
    B = [ 1/c.m, 0;
          0, 2*Caf/c.m;
          0, 2*Caf*c.lf/c.Izz;
          0, 0; 0, 0; 0, 0 ];
      
    Ad = eye(6) + A*dt;
    Bd = B*dt;
end