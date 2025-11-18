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
logged_u = out.u.Data;


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

Caf = out.Caf.Data; Car = out.Car.Data; m = 780; Izz = 1000; lf = 1.4; lr = 1.6; 
[N, ~] = size(x_smooth); % To start for loop

ddpsi_reconstructed = zeros(N,1);

for t = 1:N
    % Get the linear model (A and B) at this time step
    % We must recalculate the linear terms just as we did in the KF
    vx_safe = max(x_smooth(t,1), 1.0);
    % Yv = -(2*Caf(t) + 2*Car(t)) / (m * vx_safe(t));
    % Yr = -vx_smooth(t) - (2*Caf(t)*lf - 2*Car(t)*lr) / (m * vx_safe);
    Nv = -(2*Caf(t)*lf - 2*Car(t)*lr) / (Izz * vx_safe);
    Nr = -(2*Caf(t)*lf^2 + 2*Car(t)*lr^2) / (Izz * vx_safe);
    B_lin_psi = (2*Caf(t)*lf) / Izz;
    
    % Reconstruct the linear part of the prediction
    psi_ddot_lin = (Nv * x_smooth(t,2) + Nr * x_smooth(t,3)) + (B_lin_psi * logged_u(t,2));
    
    % Reconstruct the disturbance part of the prediction
    psi_ddot_disturbance = x_smooth(t,6) / Izz;
    
    % The "reconstructed" acceleration is the sum of both
    ddpsi_reconstructed(t) = psi_ddot_lin + psi_ddot_disturbance;
end


figure;
plot(time_vector, ddpsi_reconstructed); hold on;
plot(time_vector, out.gt_ddpsi.Data)


%% --- 4. Run the "Vehicle Ident" Optimization ---
fprintf('Data has been smoothed and is ready for optimization.\n');
fprintf('Step 1: Reconstructing ground truth target...\n');


% Prep for identification (ie optimization)
[N, ~] = size(x_smooth);
y_target = zeros(N, 1);
vx_smooth = x_smooth(:, 1);
vy_smooth = x_smooth(:, 2);
yaw_rate_smooth  = x_smooth(:, 3);
xi_psi_smooth = x_smooth(:, 6);
delta_log = logged_u(:,2); % Get the logged delta command
Caf_log = out.Caf.Data; 
Car_log = out.Car.Data; 
m = 780; Izz = 1000; lf = 1.4; lr = 1.6; 

for t = 1:N
    vx_safe = max(vx_smooth(t), 1.0);
    Caf_t = Caf_log(t);
    Car_t = Car_log(t);
    Nv = -(2*Caf_t*lf - 2*Car_t*lr) / (Izz * vx_safe);
    Nr = -(2*Caf_t*lf^2 + 2*Car_t*lr^2) / (Izz * vx_safe);
    B_lin_psi = (2*Caf_t*lf) / Izz;
    psi_ddot_lin = (Nv * vy_smooth(t) + Nr * yaw_rate_smooth(t)) + (B_lin_psi * delta_log(t));
    psi_ddot_disturbance = xi_psi_smooth(t) / Izz;
    y_target(t) = psi_ddot_lin + psi_ddot_disturbance;
end
fprintf('Ground truth target reconstructed.\n');

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
initial_guess = [20, 1.5, 1.5, 1.0, 20, 1.5, 1.5, 1.0];
lb = [0, 0.2, 0, 0, 0, 0, 0, 0];
ub = [40, 2.5, 2.5, 2.0, 40, 2.5, 2.5, 2.0];

% --- THIS IS THE FIX: Add 'Fz_data' to the function handle ---
loss_handle = @(p) loss_function(p, states_smooth, inputs_log, Fz_data, consts, y_target);
% --- END FIX ---

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxIterations', 200, 'PlotFcn', {@optimplotfval, @optimplotstepsize});

fprintf('Step 2: Starting "Vehicle Ident" optimization...\n');

[theta_optimal, final_loss] = fmincon(loss_handle, initial_guess, ...
                                      [], [], [], [], lb, ub, [], options);

% --- 5. Display Results ---
fprintf('Optimization finished.\n');
fprintf('Final Loss: %f\n', final_loss);
disp('Optimal Pacejka Parameters (theta):');
disp(theta_optimal);

disp('Ground Truth Parameters (from simulator):');
disp([22, 1.8, 1.6, 0.8, 22, 1.8, 1.6, 0.8]);
