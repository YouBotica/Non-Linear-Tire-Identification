%% --- 1. Extract the Smoothed Data for the Identification, this scrip assumes you already have x_smooth loaded in the workspace

% We need u (ie control signal of the simulation), approximate cornering
% stiffnesses and vehicle parameters
logged_u = out.u.Data;
Caf = out.Caf.Data; Car = out.Car.Data; 
m = veh_param.m; Izz = veh_param.Izz; lf = veh_param.lf; lr = veh_param.lr; 


[N, ~] = size(x_smooth); 

ddpsi_reconstructed = zeros(N,1);
ddy_reconstructed = zeros(N,1);
y_target_variance = zeros(N, 1); 

for t = 1:N
    % States & Inputs
    vx_t = x_smooth(t, 1); vy_t = x_smooth(t, 2); r_t  = x_smooth(t, 3); 
    xi_y_t = x_smooth(t, 5); xi_psi_t = x_smooth(t, 6);
    delta_t = logged_u(t, 2);
    Caf_t = out.Caf.Data(t); Car_t = out.Car.Data(t);
    
    % Linear Model Terms
    vx_safe = max(vx_t, 1.0);
    Nv = -(Caf_t*lf - Car_t*lr) / (Izz * vx_safe);
    Nr = -(Caf_t*lf^2 + Car_t*lr^2) / (Izz * vx_safe);
    Yv = -(Caf_t + Car_t) / (m * vx_safe);
    Yr = -vx_t - (Caf_t*lf - Car_t*lr) / (m * vx_safe); % Yr includes -vx
    B_lin_psi = (Caf_t*lf) / Izz;
    
    % --- A. Reconstruct Target (Inertial Lateral Accel) ---
    % Body Frame Reconstruction:
    ddy_body = Yv * vy_t + Yr * r_t + (Caf_t * 0 * delta_t / m) + (xi_y_t / m);
    
    % Inertial Frame Conversion (Matches IMU ground truth):
    % ay = dvy_body/dt + vx*r
    ddy_reconstructed(t) = ddy_body + (vx_t * r_t);

    % --- B. Reconstruct Variance ---
    P_t = P_smooth(:, :, t); % Assuming r_smooth is the smoothed covariance
    % H matrix maps state to ay = [0, Yv, (Yr+vx), 0, 1/m, 0]
    H_t = [0, Yv, (Yr + vx_t), 0, 1/m, 0];
    sigma2 = H_t * P_t * H_t';
    y_target_variance(t) = max(sigma2, 1e-9); 
end

figure; 
plot(out.x_t_t.Time, ddy_reconstructed, 'b'); hold on;
plot(out.gt_ay.Time, out.gt_ay.Data, 'r--');
legend('Reconstructed (ay)', 'Ground Truth (ay)');
title('Lateral Acceleration Check');

%% ---  Run the "Vehicle Ident" Optimization ---
fprintf('Data has been smoothed and is ready for optimization.\n');

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

% --- 1. Define Initial Guess & Bounds ---
initial_guess = [15, 1.5, 1.5, 1.0, 15, 1.5, 1.5, 1.0, 0, 0, 0];
lb = [1, 0.5, 0.1, 0.2, 1, 0.5, 0.1, 0.2, 0, 0, 0];
ub = [40, 3.0, 3.0, 2.0, 40, 3.0, 3.0, 2.0, 0, 0, 0];

states_smooth = x_smooth;
Fz_data = [out.FzF.Data, out.FzR.Data];

% --- 4. PRE-OPTIMIZATION CHECK ---
fprintf('Calculating response with Initial Guess...\n');
y_pred_initial = predict_NL_accel(initial_guess, states_smooth, logged_u, Fz_data, consts);

% --- 5. Run Optimization ---
fprintf('Optimizing against Ground Truth ay...\n');

% loss function signature: loss_function_weighted(params, states_smooth, inputs_log, Fz_data, consts, y_target, y_target_variance, initial_guess)
loss_handle = @(p) loss_function_weighted(p, states_smooth, logged_u, Fz_data, consts, ddy_reconstructed, y_target_variance, initial_guess);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxIterations', 500, 'MaxFunctionEvaluations', 10e3);

[theta_test, loss_test] = fmincon(loss_handle, initial_guess, [], [], [], [], lb, ub, [], options);

% --- 6. Results & Validation ---
disp('--- UNIT TEST RESULTS ---');
disp('Ground Truth: [22, 1.8, 1.6, 0.8]');
disp('Identified:');
disp(theta_test(1:8));

% 6a. Calculate Final Prediction (Identified Model)
y_pred_final = predict_NL_accel(theta_test, states_smooth, logged_u, Fz_data, consts);

% 6b. Calculate Ground Truth Prediction (Perfect Model Check)
% This tests if your predict_NL_accel function can reproduce the truth 
% given the correct parameters. If this doesn't match y_target_GT, 
% then predict_NL_accel has a bug.
theta_ground_truth = [22, 1.8, 1.6, 0.8, 22, 1.8, 1.6, 0.8, 0, 0, 0];
y_pred_true_params = predict_NL_accel(theta_ground_truth, states_smooth, logged_u, Fz_data, consts);

% --- 7. Combined Plot ---
figure('Name', 'Identification Result');
plot(out.gt_ay.Data(), 'k', 'LineWidth', 1.5); hold on;       % Measured Truth
plot(y_pred_initial, 'g--', 'LineWidth', 1);             % Start
plot(y_pred_final, 'r--', 'LineWidth', 1.5);             % End (Identified)
plot(y_pred_true_params, 'c:', 'LineWidth', 1.5);          % Ideal (True Params)

legend('Meas. Ground Truth (ay)', 'Initial Guess', 'Identified Model', 'True Model Prediction');
title('Unit Test Result: Optimization Progress');
grid on;

