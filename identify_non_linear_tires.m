%% Post-Processing:

% 1. Load data from Simulink
% You might need to access them like: logged_x_t_t = logsout.get('logged_x_t_t').Values;
% For this example, we assume they are already in the workspace.

% Access the .Data property from the timeseries object.
% The 'squeeze' function removes any extra (1x1xN) dimensions.
logged_x_t_t   = squeeze(out.x_t_t.Data);
logged_x_t_t_1 = squeeze(out.x_t_t_1.Data);
logged_r_t_t   = out.sigma_t_t.Data;
logged_r_t_t_1 = out.sigma_t_t_1.Data;
logged_Ad      = out.Ad.Data;


% 2. Run the RTS Smoother
[x_smooth, r_smooth] = RTS_smoother(logged_x_t_t, logged_r_t_t, ...
                                        logged_x_t_t_1, logged_r_t_t_1, ...
                                        logged_Ad);


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

Caf = 63000; Car = 63000; m = 780; Izz = 1000; lf = 1.4; lr = 1.6; 

% vx_safe = max(vx_smooth(t), 1.0);
Yv = -(2*Caf + 2*Car) / (m * x_smooth(:,1));
Yr = -x_smooth(:,1) - (2*Caf*lf - 2*Car*lr) / (m * x_smooth(:,1));
Nv = -(2*Caf*lf - 2*Car*lr) / (Izz * x_smooth(:,1));
Nr = -(2*Caf*lf^2 + 2*Car*lr^2) / (Izz * x_smooth(:,1));
B_lin_psi = (2*Caf*lf) / Izz;

ddpsi_reconstructed = T*(Nv*x_smooth(:,2) + Nr*x_smooth(:,3) + (2*Caf*lf/Izz)*out.u.Data(:,2) + x_smooth(:,6)/Izz);% + 

figure;
plot(time_vector, ddpsi_reconstructed); hold on;
plot(time_vector, out.gt_ddpsi.Data)


%%
% --- 4. Run the "Physics Whiz" Optimization ---
% ... (Optimize using fmincon or lsqnonlin
%      using the 'xi_psi_smoothed' data as the target) ...

fprintf('Data has been smoothed and is ready for optimization.\n');

fprintf('Step 1: Reconstructing ground truth target...\n');

[N, ~] = size(x_smooth);
y_target = zeros(N, 1);
vx_smooth = x_smooth(:, 1);
vy_smooth = x_smooth(:, 2);
yaw_rate_smooth  = x_smooth(:, 3);
xi_psi_smooth = x_smooth(:, 6);
delta_log = out.u.Data(:,2); % Get the logged delta command

Caf = 63000; Car = 63000; m = 780; Izz = 1000; lf = 1.4; lr = 1.6; 

for t = 1:N
    % Get the linear model (A and B) at this time step
    % We must recalculate the linear terms just as we did in the KF
    vx_safe = max(vx_smooth(t), 1.0);
    Yv = -(2*Caf + 2*Car) / (m * vx_safe);
    Yr = -vx_smooth(t) - (2*Caf*lf - 2*Car*lr) / (m * vx_safe);
    Nv = -(2*Caf*lf - 2*Car*lr) / (Izz * vx_safe);
    Nr = -(2*Caf*lf^2 + 2*Car*lr^2) / (Izz * vx_safe);
    B_lin_psi = (2*Caf*lf) / Izz;
    
    % Reconstruct the linear part of the prediction
    psi_ddot_lin = (Nv * vy_smooth(t) + Nr * yaw_rate_smooth(t)) + (B_lin_psi * delta_log(t));
    
    % Reconstruct the disturbance part of the prediction
    psi_ddot_disturbance = xi_psi_smooth(t) / Izz;
    
    % The "true" acceleration is the sum of both
    y_target(t) = psi_ddot_lin + psi_ddot_disturbance;
end


% --- 2. Prepare Data for Optimizer ---
% (We already created y_target in Step 1)
states_smooth = x_smooth;
inputs_log = out.u.Data(:,:); % Your logged [gp, delta]
consts.lf = lf;
consts.lr = lr;
consts.Izz = Izz;
consts.m = m;
consts.g = 9.81;      % Gravitational acceleration
consts.h = 0.35;       % Estimated CG height (in meters)
consts.T_w = 1.5;     % Estimated track width (in meters)
% (Add any other constants 'predict_NL_accel' needs)

% --- 3. Set up and Run Optimization ---
% Initial guess for theta = [Bf, Cf, Df, Ef, Br, Cr, Dr, Er]
% gt: B = 22; C = 1.8; D=1.6; E = 0.8;
initial_guess = [20, 1.0, 1.2, 1.5, 30, 1.9, 0.7, 1.5];

% Define lower and upper bounds 
% B(Stiffness), C(Shape), D(Peak), E(Curvature)
lb = [ 1, 1.0, 0.1, 0.1,  1, 1.0, 0.1, 0.1];
ub = [40, 2.5, 2.5, 2, 40, 2.5, 2.5, 2];

% Create an anonymous function handle for the loss
loss_handle = @(p) loss_function(p, states_smooth, inputs_log, consts, y_target);

% Set optimizer options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'MaxIterations', 100);
options.FiniteDifferenceType = 'central';
options.FiniteDifferenceStepSize = 1e-4;
options.OptimalityTolerance = 1e-6;
options.StepTolerance = 1e-8;
options.ConstraintTolerance = 1e-6;
options.PlotFcn = {@optimplotx, @optimplotfval, @optimplotstepsize};


% options = optimoptions('fmincon','Algorithm','interior-point');


fprintf('Step 2: Starting "Physics Whiz" optimization...\n');

% Run the optimizer!
[theta_optimal, final_loss] = fmincon(loss_handle, initial_guess, ...
                                      [], [], [], [], lb, ub, [], options);

% --- 4. Display Results ---
fprintf('Optimization finished.\n');
fprintf('Final Loss: %f\n', final_loss);
disp('Optimal Pacejka Parameters (theta):');
disp(theta_optimal);
