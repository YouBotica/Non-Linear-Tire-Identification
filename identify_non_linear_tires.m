% Post-Processing:

% 1. Load data from Simulink
% You might need to access them like: logged_x_t_t = logsout.get('logged_x_t_t').Values;
% For this example, we assume they are already in the workspace.

% 2. Run the RTS Smoother
[x_smooth, r_smooth] = RTS_smoother(out.x_t_t, out.sigma_t_t, ...
                                        out.x_t_t_1, out.sigma_t_t_1, ...
                                        out.Ad);

% 3. Extract the "Golden" Data
% x_smooth is an (N x 6) matrix. Let's get the smoothed disturbance:
% x_hat = [vx, vy, r, xi_x, xi_y, xi_psi]'
xi_vx_smoothed = x_smooth(:, 4);
xi_psi_smoothed = x_smooth(:, 6);
xi_vy_smoothed = x_smooth(:, 5);
vx_smoothed = x_smooth(:, 1);
vy_smoothed = x_smooth(:, 2);
yaw_rate_smoothed = x_smooth(:, 3);

% --- 4. Run the "Physics Whiz" Optimization ---
% ... (This is where you'd call your fmincon or lsqnonlin
%      using the 'xi_psi_smoothed' data as the target) ...

fprintf('Data has been smoothed and is ready for optimization.\n');
