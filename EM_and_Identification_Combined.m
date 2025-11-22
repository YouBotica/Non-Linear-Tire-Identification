%% EM_and_Identification_Combined.m
% Combined script for EM-based Kalman filter tuning and non-linear tire identification.
% This integrates the EM algorithm (from EM_EKF_tuner.m) with reconstruction and optimization (from identify_non_linear_tires.m).
% Assumes 'out' workspace variable from Simulink is loaded.

%% 1. Load and Prepare Data
if ~exist('out', 'var')
    error('Please run the Simulink simulation first to generate "out".');
end

% Extract measurements and inputs
y_log = [out.measurements.Data(:,1), out.measurements.Data(:,2), out.measurements.Data(:,3), out.measurements.Data(:,4), out.FzF.Data, out.FzR.Data];
u_log = out.u.Data;

% Vehicle constants
consts.m = veh_param.m;
consts.Izz = veh_param.Izz;
consts.lf = veh_param.lf;
consts.lr = veh_param.lr;
consts.Cd = veh_param.Cd;
consts.K_pacejka = veh_param.Bf * veh_param.Cf * veh_param.Df;  % Linear slope for dynamic stiffness
dt = control_param.T;  % Sample time

[T_steps, ~] = size(y_log);

%% 2. EM Algorithm for Q and R Tuning
fprintf('Starting EM Algorithm for Kalman Tuning...\n');

% Initial Q and R
Q_diag = control_param.ekf_Q_diag;
Q_curr = diag(Q_diag);
R_curr = diag([control_param.ekf_vx_sensor_noise; control_param.ekf_vy_sensor_noise; control_param.ekf_dpsi_sensor_noise; control_param.ekf_ay_sensor_noise]);

MAX_ITER = 12;  % Typical convergence
for iter = 1:MAX_ITER
    fprintf('  EM Iteration %d...\n', iter);
    
    % E-Step: Forward EKF + RTS Smoother
    [x_filt, P_filt, x_pred, P_pred, Ad_log] = forward_KF_EM(y_log, u_log, Q_curr, R_curr, control_param.ekf_initial_cov, consts, dt);
    [x_smooth, P_smooth, P_cross_smooth] = RTS_smoother_EM(x_filt, P_filt, x_pred, P_pred, Ad_log);
    
    % M-Step: Update Q and R
    % Update R
    R_new_sum = zeros(4, 4);
    for t = 1:T_steps
        xt = x_smooth(t, :)';
        Pt = P_smooth(:, :, t);
        yt = y_log(t, 1:4)';
        Fz = [out.FzF.Data(t); out.FzR.Data(t)];
        [Ad, Bd, Yv, Yr] = get_linear_model(xt, Fz, consts, dt);
        C_aug = [eye(3), zeros(3); 0, Yv, (Yr + xt(1)), 0, 1/consts.m, 0];
        res = yt - C_aug * xt;
        R_new_sum = R_new_sum + (res * res' + C_aug * Pt * C_aug');
    end
    R_calculated = R_new_sum / T_steps;
    
    % Update Q
    Q_new_sum = zeros(6, 6);
    for t = 2:T_steps
        xt = x_smooth(t, :)';
        xt_prev = x_smooth(t-1, :)';
        ut_prev = u_log(t-1, :)';
        Pt = P_smooth(:, :, t);
        Pt_prev = P_smooth(:, :, t-1);
        P_cross = P_cross_smooth(:, :, t);
        [Ad, Bd, Yv, Yr] = get_linear_model(xt_prev, y_log(t-1,5:6), consts, dt);
        pred_res = xt - (Ad * xt_prev + Bd * ut_prev);
        cov_correction = Pt + (Ad * Pt_prev * Ad') - (P_cross * Ad') - (Ad * P_cross');
        Q_new_sum = Q_new_sum + (pred_res * pred_res' + cov_correction);
    end
    Q_calculated = Q_new_sum / (T_steps - 1);
    
    % Apply constraints
    Q_new = diag(diag(Q_calculated));
    R_new = diag(diag(R_calculated));
    Q_new(1:3,1:3) = diag(Q_diag(1:3));  % Lock physics states
    SAFE_FLOOR = 1e-9;
    Q_new = max(Q_new, eye(6) * SAFE_FLOOR);
    R_new = max(R_new, eye(4) * SAFE_FLOOR);
    
    Q_curr = Q_new;
    R_curr = R_new;  % Enable R updates for better EM convergence

    % Compute Log-Likelihood for Monitoring
    log_likelihood = 0;
    for t = 1:T_steps
        xt = x_smooth(t, :)';
        Pt = P_smooth(:, :, t);
        yt = y_log(t, 1:4)';
        Fz = [out.FzF.Data(t); out.FzR.Data(t)];
        [Ad, Bd, Yv, Yr] = get_linear_model(xt, Fz, consts, dt);
        C_aug = [eye(3), zeros(3); 0, Yv, (Yr + xt(1)), 0, 1/consts.m, 0];
        innov = yt - C_aug * xt;
        innov_cov = C_aug * Pt * C_aug' + R_curr;
        log_likelihood = log_likelihood - 0.5 * (innov' / innov_cov * innov + log(det(innov_cov)) + 4*log(2*pi));
    end
    log_likelihood_history(iter) = log_likelihood;
    fprintf('    Log-Likelihood: %.4f\n', log_likelihood);
end

fprintf('EM Converged. Learned Q: %s\n', mat2str(diag(Q_curr)));
fprintf('Learned R: %s\n', mat2str(diag(R_curr)));
fprintf('Log-Likelihood History: %s\n', mat2str(log_likelihood_history));

fprintf('EM Converged. Learned Q: %s\n', mat2str(diag(Q_curr)));
fprintf('Learned R: %s\n', mat2str(diag(R_curr)));

%% 3. Reconstruction and Optimization
fprintf('Starting Reconstruction and Optimization...\n');

% Reconstruct lateral acceleration using EM-tuned smoothed states
ddy_reconstructed = zeros(T_steps, 1);
y_target_variance = zeros(T_steps, 1);

for t = 1:T_steps
    vx_t = x_smooth(t, 1); vy_t = x_smooth(t, 2); r_t = x_smooth(t, 3);
    xi_y_t = x_smooth(t, 5);
    delta_t = u_log(t, 2);
    FzF_t = y_log(t, 5); FzR_t = y_log(t, 6);
    Caf_t = consts.K_pacejka * FzF_t;
    Car_t = consts.K_pacejka * FzR_t;
    
    vx_safe = max(vx_t, 1.0);
    Yv = -(Caf_t + Car_t) / (consts.m * vx_safe);
    Yr = -vx_t - (Caf_t*consts.lf - Car_t*consts.lr) / (consts.m * vx_safe);
    
    ddy_body = Yv * vy_t + Yr * r_t + (Caf_t * delta_t / consts.m) + (xi_y_t / consts.m);
    ddy_reconstructed(t) = ddy_body + vx_t * r_t;  % Inertial frame
    
    % Variance
    P_t = P_smooth(:, :, t);
    H_t = [0, Yv, Yr, 0, 1/consts.m, 0];
    sigma2 = H_t * P_t * H_t';
    y_target_variance(t) = max(sigma2, 1e-9);
end

% Optimization setup
initial_guess = [20, 1.5, 1.5, 1.0, 20, 1.5, 1.5, 1.0, 0, 0, 0];
lb = [1, 0.5, 0.1, 0.2, 1, 0.5, 0.1, 0.2, -100, -100, -100];
ub = [40, 3.0, 3.0, 2.0, 40, 3.0, 3.0, 2.0, 100, 100, 100];

problem = createOptimProblem('fmincon', ...
    'objective', @(p) loss_function_weighted(p, x_smooth, u_log, [y_log(:,5), y_log(:,6)], consts, ddy_reconstructed, y_target_variance, initial_guess), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 2000));

ms = MultiStart('Display', 'final', 'UseParallel', false);
[theta_optimal, final_loss] = run(ms, problem, 5);

fprintf('Optimization Complete. Final Loss: %.4f\n', final_loss);
disp('Optimal Pacejka Params:');
disp(theta_optimal(1:8));
disp('Optimal Residuals:');
disp(theta_optimal(9:11));

%% Helper Functions (Inline for Self-Containment)

function [x_hist, P_hist, x_pred_hist, P_pred_hist, Ad_hist] = forward_KF_EM(y, u, Q, R, P_initial, c, dt)
    [N, ~] = size(y);
    x = zeros(6, 1); x(1) = 5.0; P = P_initial;
    x_hist = zeros(N, 6); P_hist = zeros(6,6,N); x_pred_hist = zeros(N,6); P_pred_hist = zeros(6,6,N); Ad_hist = zeros(6,6,N);
    for k = 1:N
        uk = u(k, :)'; yk = y(k, 1:4)'; Fz = y(k, 5:6);
        [Ad, Bd, Yv, Yr] = get_linear_model(x, Fz, c, dt);
        C_aug = [eye(3), zeros(3); 0, Yv, (Yr + x(1)), 0, 1/c.m, 0];
        x_pred = Ad * x + Bd * uk; P_pred = Ad * P * Ad' + Q;
        K = P_pred * C_aug' / (C_aug * P_pred * C_aug' + R);
        x = x_pred + K * (yk - C_aug * x_pred);
        P = (eye(6) - K*C_aug) * P_pred * (eye(6) - K*C_aug)' + K*R*K';
        x_hist(k,:) = x'; P_hist(:,:,k) = P; x_pred_hist(k,:) = x_pred'; P_pred_hist(:,:,k) = P_pred; Ad_hist(:,:,k) = Ad;
    end
end

function [x_smooth, P_smooth, P_cross_smooth] = RTS_smoother_EM(x_filt, P_filt, x_pred, P_pred, Ad_hist)
    [N, num_states] = size(x_filt);
    x_smooth = zeros(N, num_states); P_smooth = zeros(num_states, num_states, N); P_cross_smooth = zeros(num_states, num_states, N);
    x_smooth(N,:) = x_filt(N,:); P_smooth(:,:,N) = P_filt(:,:,N);
    for k = (N-1):-1:1
        x_k_k = x_filt(k,:)'; P_k_k = P_filt(:,:,k);
        x_k1_k = x_pred(k+1,:)'; P_k1_k = P_pred(:,:,k+1); Ad_k = Ad_hist(:,:,k);
        jitter = 1e-9 * eye(size(P_k1_k));
        J = (P_k_k * Ad_k') / (P_k1_k + jitter);
        x_smooth(k,:) = (x_k_k + J * (x_smooth(k+1,:)' - x_k1_k))';
        P_smooth(:,:,k) = P_k_k + J * (P_smooth(:,:,k+1) - P_k1_k) * J';
        P_cross_smooth(:,:,k+1) = P_smooth(:,:,k+1) * J';
    end
end

function [Ad, Bd, Yv, Yr] = get_linear_model(x, Fz, c, dt)
    vx = max(x(1), 1.0); FzF = Fz(1); FzR = Fz(2);
    Caf = c.K_pacejka * FzF; Car = c.K_pacejka * FzR;
    Yv = -(Caf + Car) / (c.m * vx); Yr = -vx - (Caf*c.lf - Car*c.lr) / (c.m * vx);
    Nv = -(Caf*c.lf - Car*c.lr) / (c.Izz * vx); Nr = -(Caf*c.lf^2 + Car*c.lr^2) / (c.Izz * vx);
    A = [-c.Cd/c.m, 0, 0, 1/c.m, 0, 0; 0, Yv, Yr, 0, 1/c.m, 0; 0, Nv, Nr, 0, 0, 1/c.Izz; zeros(3,6)];
    B = [1/c.m, 0; 0, Caf/c.m; 0, Caf*c.lf/c.Izz; zeros(3,2)];
    Ad = eye(6) + A*dt; Bd = B*dt;
end

function loss = loss_function_weighted(params, states, inputs, Fz_data, consts, y_target, y_target_variance, initial_guess)
    y_predict = predict_NL_accel(params, states, inputs, Fz_data, consts);
    residuals = y_target - y_predict;
    safe_variance = max(y_target_variance, 1e-6);
    data_loss = sum((residuals.^2) ./ safe_variance);
    lambda_prior = 0.1;
    scale_factors = abs(initial_guess); scale_factors(scale_factors < 1e-6) = 1.0;
    prior_weights = [1,1,1,1,1,1,1,1,0.1,0.1,0.1];
    param_dev = ((params - initial_guess) ./ scale_factors) .* prior_weights;
    prior_loss = lambda_prior * sum(param_dev.^2);
    lambda_sim = 0.1;  % Reduced to allow asymmetry
    p_front = params(1:4); p_rear = params(5:8);
    sim_scale = abs(initial_guess(1:4));
    sim_diff = (p_front - p_rear) ./ sim_scale;
    similarity_loss = lambda_sim * sum(sim_diff.^2);
    loss = data_loss + prior_loss + similarity_loss;
    if isnan(loss) || isinf(loss), loss = 1e20; end
end

function y_ddot_nl = predict_NL_accel(params, states, inputs, Fz_data, consts)
    Bf = params(1); Cf = params(2); Df = params(3); Ef = params(4);
    Br = params(5); Cr = params(6); Dr = params(7); Er = params(8);
    lf = consts.lf; lr = consts.lr; Izz = consts.Izz; m = consts.m;
    vx = states(:,1); vy = states(:,2); r = states(:,3);
    delta = inputs(:,2);
    Fz_f_dyn = Fz_data(:,1); Fz_r_dyn = Fz_data(:,2);
    vx_safe = max(vx, 1.0);
    alpha_f = atan((vy + lf.*r)./vx_safe) - delta;
    alpha_r = atan((vy - lr.*r)./vx_safe);
    F_peak_f = Df .* Fz_f_dyn; x_f = Bf .* alpha_f;
    Fyf = -F_peak_f .* sin(Cf .* atan(x_f - Ef .* (x_f - atan(x_f))));
    F_peak_r = Dr .* Fz_r_dyn; x_r = Br .* alpha_r;
    Fyr = -F_peak_r .* sin(Cr .* atan(x_r - Er .* (x_r - atan(x_r))));
    y_ddot_nl = -vx.*r + (Fyf + Fyr) ./ m;
end