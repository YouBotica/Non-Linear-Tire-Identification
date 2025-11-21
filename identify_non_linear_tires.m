%% Post-Processing:

% 1. Load data from Simulink
logged_x_t_t   = squeeze(out.x_t_t.Data);
logged_x_t_t_1 = squeeze(out.x_t_t_1.Data);
logged_r_t_t   = out.sigma_t_t.Data;
logged_r_t_t_1 = out.sigma_t_t_1.Data;
logged_Ad      = out.Ad.Data;
logged_u = out.u.Data;


% % 2. Run the RTS Smoother
% [x_smooth, r_smooth] = RTS_smoother(logged_x_t_t, logged_r_t_t, ...
%                                         logged_x_t_t_1, logged_r_t_t_1, ...
%                                         logged_Ad); DEPRECATED as we run
%                                         EM to get optimal params.
%                                        


% 3. Extract the "Golden" Data
% x_smooth is an (N x 6) matrix. Let's get the smoothed disturbance:
% x_hat = [vx, vy, r, xi_x, xi_y, xi_psi]'
xi_vx_smoothed = x_smooth(:, 4);
xi_psi_smoothed = x_smooth(:, 6);
xi_vy_smoothed = x_smooth(:, 5);
vx_smoothed = x_smooth(:, 1);
vy_smoothed = x_smooth(:, 2);
yaw_rate_smoothed = x_smooth(:, 3);

% --- Plot Smoother vs. Filter Results ---
% This script assumes 'logged_x_t_t' (filtered) and 'x_smooth' (smoothed)
% are in your MATLAB workspace.

fprintf('Plotting Filter vs. Smoother states...\n');

% --- 1. Define Constants ---
[N, num_states] = size(logged_x_t_t);
T = 0.01; % Your simulation sample time
time_vector = (0:N-1) * T;

% Define the names for your states
state_names = {
    'vx (m/s)', ...
    'vy (m/s)', ...
    'r (rad/s)', ...
    'xi_x (Disturbance)', ...
    'xi_y (Disturbance)', ...
    'xi_psi (Disturbance)'
};

% --- 2. Create the Figure ---
figure;
sgtitle('Kalman Filter (Real-Time) vs. RTS Smoother (Offline)', 'FontSize', 14);

% --- 3. Loop and Plot ---
for i = 1:num_states
    % Create a 3x2 grid of subplots
    subplot(3, 2, i);
    
    % Plot the real-time KF estimate (orange, thin)
    plot(time_vector, logged_x_t_t(:, i), 'r-', 'LineWidth', 1);
    
    hold on;
    
    % Plot the offline Smoother estimate (blue, thick)
    plot(time_vector, x_smooth(:, i), 'b-', 'LineWidth', 2);
    
    % Make it look good
    title(state_names{i});
    xlabel('Time (s)');
    ylabel('Value');
    grid on;
    legend('Filter (x_t_t)', 'Smoother (x_t|T)');
    hold off;
end

fprintf('Plotting complete.\n');

% Test estimated angular acceleration vs ground truth from Simulink's
% Vehicle Dynamics Blockset:

Caf = out.Caf.Data; Car = out.Car.Data; m = 780; Izz = 1000; lf = 1.4; lr = 1.6; 
[N, ~] = size(x_smooth); % To start for loop

ddpsi_reconstructed = zeros(N,1);
y_target_variance = zeros(N, 1); % Pre-allocate variance

for t = 1:N
    % States & Inputs
    vx_t = x_smooth(t, 1); vy_t = x_smooth(t, 2); r_t  = x_smooth(t, 3); xi_psi_t = x_smooth(t, 6);
    delta_t = logged_u(t, 2);
    Caf_t = out.Caf.Data(t); Car_t = out.Caf.Data(t);
    
    % Linear Model Terms
    vx_safe = max(vx_t, 1.0);
    Nv = -(2*Caf_t*lf - 2*Car_t*lr) / (Izz * vx_safe);
    Nr = -(2*Caf_t*lf^2 + 2*Car_t*lr^2) / (Izz * vx_safe);
    B_lin_psi = (2*Caf_t*lf) / Izz;
    
    % A. Reconstruct Target (Mean)
    psi_ddot_lin = (Nv * vy_t + Nr * r_t) + (B_lin_psi * delta_t);
    psi_ddot_dist = xi_psi_t / Izz;
    ddpsi_reconstructed(t) = psi_ddot_lin + psi_ddot_dist;
    
    % B. Reconstruct Variance (Uncertainty)
    % Var(y) = H * P * H'
    P_t = P_smooth(:, :, t);
    H_t = [0, Nv, Nr, 0, 0, 1/Izz]; % Linear mapping vector
    
    sigma2 = H_t * P_t * H_t';
    y_target_variance(t) = max(sigma2, 1e-9); % Floor
end


figure;
plot(time_vector, ddpsi_reconstructed); hold on;
plot(time_vector, out.gt_ddpsi.Data)


%% --- 4. Run the "Vehicle Ident" Optimization ---
fprintf('Data has been smoothed and is ready for optimization.\n');
fprintf('Step 1: Reconstructing ground truth target...\n');


% --- 2. Prepare Data for Optimizer ---
states_smooth = x_smooth;
inputs_log = logged_u; % Your logged [gp, delta]

% --- THIS IS THE FIX: Load Fz data as a separate variable ---
% (Make sure your 'To Workspace' blocks are named 'FzF_log' and 'FzR_log')
Fz_data = [out.FzF.Data, out.FzR.Data]; % Creates an N x 2 matrix
% --- END FIX ---

% --- 3. Populate consts struct ---
consts.lf = lf;
consts.lr = lr;
consts.Izz = Izz;
consts.m = m;
consts.g = 9.81;      % Gravitational acceleration
consts.h = 0.35;      % CG height (Make sure this matches your Simulink block!)
% (We removed T_w as it's not in your sim's equations)

% --- 4. Set up and Run Optimization ---
% [Bf, Cf, Df, Ef, Br, Cr, Dr, Er, P1, P2, P3]
initial_guess = [10, 1.5, 1.5, 1.0, 1, 2.0, 0.9, 1.0, 0, 0, 0];

% Bounds (Pacejka > 0, Residuals unbounded)
lb = [1, 0.5, 0.1, 0.2, 1, 0.5, 0.1, 0.2, -100, -100, -100];
ub = [40, 3.0, 3.0, 2.0, 40, 3.0, 3.0, 2.0,  100, 100, 100];

% Create Function Handle (Passes all data + variance + initial_guess)
loss_handle = @(p) loss_function_weighted(p, states_smooth, inputs_log, Fz_data, consts, ddpsi_reconstructed, y_target_variance, initial_guess);

% --- 5. Run MultiStart ---
fprintf('Step 3: Running MultiStart...\n');
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
                       'MaxIterations', 5000, 'StepTolerance', 1e-8);

problem = createOptimProblem('fmincon', ...
    'objective',   loss_handle, ...
    'x0',          initial_guess, ...
    'lb',          lb, ...
    'ub',          ub, ...
    'options',     options);

ms = MultiStart('Display', 'iter', 'UseParallel', true); % Set to false if crashing
[theta_optimal, final_loss] = run(ms, problem, 50); 

% --- 6. Results ---
fprintf('\nOptimization Finished.\n');
fprintf('Final Weighted Loss: %f\n', final_loss);

disp('--- Optimal Pacejka (1-8) ---');
disp(theta_optimal(1:8));
disp('--- Optimal Residuals (9-11) ---');
disp(theta_optimal(9:11));

disp('--- Ground Truth (Simulator) ---');
disp([22, 1.8, 1.6, 0.8, 22, 1.8, 1.6, 0.8]);