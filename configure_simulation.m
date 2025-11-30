% --- 1. Define Constants ---
%%
% Controller sampling time:
control_param.T = 0.01;

control_param.vx_noise_pwr = 0.0001;
control_param.vy_noise_pwr = 0.0000001;
control_param.dpsi_noise_pwr = 0.0000001; 

control_param.ekf_vx_sensor_noise = 0.0093; 
control_param.ekf_vy_sensor_noise = 5.8534e-5;
control_param.ekf_dpsi_sensor_noise = 6.1532e-6;
control_param.ekf_ay_sensor_noise = 0.0026;
control_param.ekf_initial_cov = [diag([1e-5; 1e-3; 1e-3; 100; 100; 100]*0)]; % Start with zero?
control_param.deccel_threshold = 1;
control_param.K1_vx = 4;
control_param.K1_dpsi = 30; 

control_param.lookahead_distance = 0.6;


control_param.ekf_Q_diag = [
    4.9439e-4;  % vx (Trust model)
    5.8044e-4;  % vy (Trust model)
    5.55501e-6;  % r (Trust model)
    17.535;  % xi_x (Disturbance can change a lot!)
    1.3216e3;  % xi_y (Disturbance can change a lot!)
    28.569; % xi_psi (Disturbance can change a lot!)
];


% Vehicle constants (to the best of our believe):
veh_param.Caf = 63; % Approximated cornering stiffness of the front tire
veh_param.Car = 63; % Approximated cornering stiffness of the rear tire
veh_param.lf = 1.4; % Distance from CG to front axle
veh_param.lr = 1.6; % Distance from CG to rear axle
veh_param.m = 780; % Vehicle's mass (780 kg)
veh_param.Izz = 1000; % Vehicle z (vertical) moment of inertia kg*m^2
veh_param.Cd = 0.3; % Vehicle Aero Coefficient
veh_param.g = 9.81;
veh_param.h = 0.35;
veh_param.initial_states = [5; 0; 0; 0; 0; 0]; % vx, control_param.dpsi_noise_pwrvy, dpsi, disturbance_x, disturbance_y, disturbance_dpsi

% Front and rear tires Magic Formula:
% These are the ground truth (unobserved) parameters of the non-linear tire
% model we are trying to estimate
veh_param.Br = 22; 
veh_param.Cr = 1.8;
veh_param.Dr = 1.6;
veh_param.Er = 0.8;
veh_param.Bf = 22; 
veh_param.Cf = 1.8;
veh_param.Df = 1.6;
veh_param.Ef = 0.8;


% --- 2. Define Track Parameters ---
longStraight = 2500;   % meters
shortStraight = 1000;  % meters
turnRadius = 600;      % meters
pointsPerSegment = 1000; % 100 points for each of the 8 segments
v_straight = 20;
v_turn = 20;

% --- 3. Generate the Track ---
fprintf('Generating roval track...\n');
track = createOvalTrack(longStraight, shortStraight, turnRadius, ...
    pointsPerSegment, v_straight, v_turn);
fprintf('Track generated successfully.\n');

% --- 4. Plot the Track ---
fprintf('Plotting track...\n');
figure;
hold on;
plot(track.X, track.Y, 'b-', 'LineWidth', 2, 'DisplayName', 'Track Path');

% Plot heading arrows
step = 50; 
quiver(track.X(1:step:end), track.Y(1:step:end), ...
       cos(track.Yaw(1:step:end)), sin(track.Yaw(1:step:end)), ...
       20, 'r'); % '20' scales the arrow length

% --- 5. Make the Plot Look Good ---
axis equal; % <-- Essential!
title('Generated "Roval" Track (Indy-style)');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Track Path', 'Track Heading');
grid on;
hold off;