%% Post-Processing:
% 1. Load data from Simulink
fprintf('Loading data from Simulink...\n');
logged_x_t_t   = squeeze(out.x_t_t.Data);
logged_x_t_t_1 = squeeze(out.x_t_t_1.Data);
logged_r_t_t   = out.sigma_t_t.Data;
logged_r_t_t_1 = out.sigma_t_t_1.Data;
logged_Ad      = out.Ad.Data;
logged_u       = out.u.Data;
logged_Caf     = out.Caf.Data;
logged_Car     = out.Car.Data;
logged_FzF     = out.FzF.Data;
logged_FzR     = out.FzR.Data;
logged_gt_ddpsi = out.gt_ddpsi.Data;

% 2. Run the RTS Smoother
[x_smooth, r_smooth] = RTS_smoother(logged_x_t_t, logged_r_t_t, ...
                                        logged_x_t_t_1, logged_r_t_t_1, ...
                                        logged_Ad);

% 3. Extract the "Golden" Data (for plotting)
fprintf('Plotting Filter vs. Smoother states...\n');
[N, num_states] = size(logged_x_t_t);
T = 0.01;
time_vector = (0:N-1) * T;
state_names = {'vx (m/s)', 'vy (m/s)', 'r (rad/s)', ...
               'xi_x (Disturbance)', 'xi_y (Disturbance)', 'xi_psi (Disturbance)'};
figure;
sgtitle('Kalman Filter (Real-Time) vs. RTS Smoother (Offline)', 'FontSize', 14);
for i = 1:num_states
    subplot(3, 2, i);
    plot(time_vector, logged_x_t_t(:, i), 'r-', 'LineWidth', 1);
    hold on;
    plot(time_vector, x_smooth(:, i), 'b-', 'LineWidth', 2);
    title(state_names{i});
    xlabel('Time (s)');
    ylabel('Value');
    grid on;
    legend('Filter (x_t_t)', 'Smoother (x_t|T)');
    hold off;
end
fprintf('Plotting complete.\n');

%% --- 4. Run the "Vehicle Ident" (Physics Whiz) Optimization ---
fprintf('Step 1: Reconstructing ground truth target (y_target)...\n');

% --- THIS IS THE *ONLY* RECONSTRUCTION LOOP ---
% We build our 'y_target' vector from the smoothed data.
% This is the "Data" from our q_phi Recognition Model.

% Define constants
m = 780; Izz = 1000; lf = 1.4; lr = 1.6; 
y_target = zeros(N,1); % This is our 'ddpsi_reconstructed'
y_target_variance = zeros(N, 1);

for t = 1:N
    % Get all smoothed states for this time step
    vx_t = x_smooth(t, 1);
    vy_t = x_smooth(t, 2);
    r_t  = x_smooth(t, 3);
    xi_psi_t = x_smooth(t, 6);
    
    % Get dynamic stiffness and input for this time step
    Caf_t = logged_Caf(t);
    Car_t = logged_Car(t);
    delta_t = logged_u(t, 2);
    
    % Recalculate the linear terms (must match the KF)
    vx_safe = max(vx_t, 1.0);
    Nv = -(2*Caf_t*lf - 2*Car_t*lr) / (Izz * vx_safe);
    Nr = -(2*Caf_t*lf^2 + 2*Car_t*lr^2) / (Izz * vx_safe);
    B_lin_psi = (2*Caf_t*lf) / Izz;
    
    % Reconstruct the linear part of the prediction
    psi_ddot_lin = (Nv * vy_t + Nr * r_t) + (B_lin_psi * delta_t);
    
    % Reconstruct the disturbance part of the prediction
    psi_ddot_disturbance = xi_psi_t / Izz;
    
    % The "true" acceleration is the sum of both
    y_target(t) = psi_ddot_lin + psi_ddot_disturbance;

    % Calculate the variance of y_target:
    % Get the 6x6 smoothed covariance matrix for this time step
    r_t_T = r_smooth(:, :, t);
    
    % Define the H-vector that maps our state [vx,vy,r,xi_x,xi_y,xi_psi]
    % to our y_target (which is Nv*vy + Nr*r + 1/Izz * xi_psi)
    H_t = [0, Nv, Nr, 0, 0, 1/Izz];
    
    % The variance of a linear combo is Var(Hx) = H * P * H'
    % This is the "variance sandwich rule"
    sigma_t_squared = H_t * r_t_T * H_t';
    
    % Store it (with a "jitter" floor to prevent division by zero)
    y_target_variance(t) = max(sigma_t_squared, 1e-9);
end

fprintf('Ground truth target reconstructed.\n');

% --- (SANITY CHECK PLOT) ---
figure;
plot(time_vector, y_target, 'b--', 'LineWidth', 1.5);
hold on;
plot(time_vector, logged_gt_ddpsi, 'r-');
legend('Reconstructed (y\_target)', 'Ground Truth (from Sim)');
title('Sanity Check: Reconstructed vs. Ground Truth');
grid on;
hold off;


% --- 5. Prepare Data for Optimizer ---
fprintf('Step 2: Preparing data for optimization...\n');
states_smooth = x_smooth; % N x 6
inputs_log = logged_u; % N x 2 ([gp, delta])
Fz_data = [logged_FzF, logged_FzR]; % N x 2 ([FzF, FzR])

consts.lf = lf;
consts.lr = lr;
consts.Izz = Izz;
consts.m = m;
consts.g = 9.81;
consts.h = 0.35; % Make sure this matches your sim block!

% --- 6. Set up and Run Optimization with MultiStart ---
fprintf('Step 3: Starting "Vehicle Ident" (Grey Box) optimization with MultiStart...\n');

% --- NEW 11-PARAMETER GUESSES & BOUNDS ---
% [Bf, Cf, Df, Ef, Br, Cr, Dr, Er, P1, P2, P3]

% Ground Truth (for first 8): [22, 1.8, 1.6, 0.8, 22, 1.8, 1.6, 0.8]
initial_guess = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];                             % Residual (start at 0)

lb = [ 5, 0, 0, 0, 10, 0, 0, 0, ... % Pacejka
      0, 0, 0 ];                           % Residual (let it be negative or positive)
      
ub = [ 40, 2.5, 2.5, 2.0, 40, 2.5, 2.5, 2.0, ... % Pacejka
       0, 0, 0 ];                           % Residual

% --- (This part is UNCHANGED, but now uses the 11-param vectors) ---
loss_handle = @(p) loss_function(p, states_smooth, inputs_log, Fz_data, consts, y_target, y_target_variance, initial_guess);

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp'); 
problem = createOptimProblem('fmincon', ...
    'objective',   loss_handle, ...
    'x0',          initial_guess*0, ...
    'lb',          lb, ...
    'ub',          ub, ...
    'options',     options);

ms = MultiStart('Display', 'iter', 'UseParallel', true);
[theta_optimal, final_loss] = run(ms, problem, 10); 

% --- 7. Display Results (Updated for 11 params) ---
fprintf('Optimization finished.\n');
fprintf('Final Loss: %f\n', final_loss);

disp('Optimal Pacejka Parameters (theta 1-8):');
disp(theta_optimal(1:8));

disp('Optimal Residual Parameters (theta 9-11):');
disp(theta_optimal(9:11));

disp('Ground Truth Parameters (from simulator):');
disp([22, 1.8, 1.6, 0.8, 22, 1.8, 1.6, 0.8]);