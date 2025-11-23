%% --- UNIT TEST: Ground Truth Identification ---
fprintf('--- STARTING GROUND TRUTH UNIT TEST ---\n');

% --- 1. Load Data & Setup Constants ---
consts.m = veh_param.m;
consts.Izz = veh_param.Izz;
consts.lf = veh_param.lf;
consts.lr = veh_param.lr;
consts.Cd = veh_param.Cd;
consts.K_pacejka = veh_param.Bf * veh_param.Cf * veh_param.Df; 
consts.g = veh_param.g;
consts.h = veh_param.h;
dt = control_param.T;

% Extract Raw Measurements (Y) and Inputs (U)
y_log = [out.measurements.Data(:,1), out.measurements.Data(:,2), ...
         out.measurements.Data(:,3), out.measurements.Data(:,4), ...
         out.FzF.Data, out.FzR.Data];
y_target_GT = y_log(:, 4); % Lateral Accel (ay) from IMU

% Inputs log:
u_log = out.u.Data;

% EKF/Smoother Setup (Standard)
Q_diag = control_param.ekf_Q_diag;
Q_curr = diag(Q_diag);
R_curr = diag([control_param.ekf_vx_sensor_noise; ...
               control_param.ekf_vy_sensor_noise; ...
               control_param.ekf_dpsi_sensor_noise; ...
               control_param.ekf_ay_sensor_noise]);
P_initial = control_param.ekf_initial_cov;

% Run Forward Pass
[x_human_filt, P_human_filt, x_human_pred, P_human_pred, Ad_log] = ...
    forward_KF_EM(y_log, u_log, Q_curr, R_curr, P_initial, consts, dt);

% Run Backward Pass (Smoother)
[x_smooth, ~, ~] = RTS_smoother_EM(x_human_filt, P_human_filt, x_human_pred, P_human_pred, Ad_log);

% --- 2. Setup Optimization Data ---
N = length(y_target_GT);
y_target_var_GT = ones(N, 1) * 1e-4; 
states_smooth = x_smooth;
inputs_log = u_log;
Fz_data = [out.FzF.Data, out.FzR.Data]; 

% --- 3. Define Initial Guess & Bounds ---
initial_guess = [15, 1.5, 1.5, 1.0, 15, 1.5, 1.5, 1.0, 0, 0, 0];
lb = [1, 0.5, 0.1, 0.2, 1, 0.5, 0.1, 0.2, 0, 0, 0];
ub = [40, 3.0, 3.0, 2.0, 40, 3.0, 3.0, 2.0, 0, 0, 0];

% --- 4. PRE-OPTIMIZATION CHECK ---
fprintf('Calculating response with Initial Guess...\n');
y_pred_initial = predict_NL_accel(initial_guess, states_smooth, inputs_log, Fz_data, consts);

% --- 5. Run Optimization ---
fprintf('Optimizing against Ground Truth ay...\n');
loss_handle = @(p) loss_function_weighted(p, states_smooth, inputs_log, Fz_data, consts, y_target_GT, y_target_var_GT, initial_guess);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxIterations', 500, 'MaxFunctionEvaluations', 10e3);

[theta_test, loss_test] = fmincon(loss_handle, initial_guess, [], [], [], [], lb, ub, [], options);

% --- 6. Results & Validation ---
disp('--- UNIT TEST RESULTS ---');
disp('Ground Truth: [22, 1.8, 1.6, 0.8]');
disp('Identified:');
disp(theta_test(1:8));

% 6a. Calculate Final Prediction (Identified Model)
y_pred_final = predict_NL_accel(theta_test, states_smooth, inputs_log, Fz_data, consts);

% 6b. Calculate Ground Truth Prediction (Perfect Model Check)
% This tests if your predict_NL_accel function can reproduce the truth 
% given the correct parameters. If this doesn't match y_target_GT, 
% then predict_NL_accel has a bug.
theta_ground_truth = [22, 1.8, 1.6, 0.8, 22, 1.8, 1.6, 0.8, 0, 0, 0];
y_pred_true_params = predict_NL_accel(theta_ground_truth, states_smooth, inputs_log, Fz_data, consts);

% --- 7. Combined Plot ---
figure('Name', 'Identification Result');
plot(y_target_GT, 'k', 'LineWidth', 2.5); hold on;       % Measured Truth
plot(y_pred_initial, 'g--', 'LineWidth', 1);             % Start
plot(y_pred_final, 'r--', 'LineWidth', 1.5);             % End (Identified)
plot(y_pred_true_params, 'c:', 'LineWidth', 2);          % Ideal (True Params)

legend('Meas. Ground Truth (ay)', 'Initial Guess', 'Identified Model', 'True Model Prediction');
title('Unit Test Result: Optimization Progress');
grid on;