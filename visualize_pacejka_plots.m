%% Compare Identified Tire Models vs. Ground Truth
clc; close all;

% --- 1. Define Parameters [B, C, D, E] ---
theta_gt    = [11,      2.0,    1.0,    0.7];
theta_front = [17.3461, 1.2192, 2.6417, 0.7022];
theta_rear  = [24.4151, 1.7344, 1.7046, 1.4463];

% --- 2. Simulation Constants ---
m = 780;     % Mass (kg)
g = 9.81;    % Gravity (m/s^2)
Fz_nominal = 40000; %(m * g) / 2; % Nominal load per axle/tire (~3826 N)

% --- 3. Generate Slip Angle Vector ---
% Plot from -0.3 to 0.3 radians (approx -17 to +17 degrees)
alpha = linspace(-0.3, 0.3, 1000);

% --- 4. Define Pacejka Function ---
% Fy = (D*Fz) * sin(C * atan(B*alpha - E*(B*alpha - atan(B*alpha))))
pacejka_fun = @(p, a) (p(3) * Fz_nominal) .* ...
    sin(p(2) .* atan(p(1).*a - p(4).*(p(1).*a - atan(p(1).*a))));

% --- 5. Calculate Forces ---
Fy_gt    = pacejka_fun(theta_gt, alpha);
Fy_front = pacejka_fun(theta_front, alpha);
Fy_rear  = pacejka_fun(theta_rear, alpha);

% --- 6. Plotting ---
figure('Name', 'Tire Identification Check', 'Color', 'w', 'Position', [100 100 800 600]);

plot(alpha, Fy_gt, 'k', 'LineWidth', 3, 'DisplayName', 'Ground Truth'); 
hold on;
plot(alpha, Fy_front, 'r--', 'LineWidth', 2, 'DisplayName', 'Identified Front');
plot(alpha, Fy_rear, 'b-.', 'LineWidth', 2, 'DisplayName', 'Identified Rear');

% Formatting
grid on;
xlabel('Slip Angle \alpha (rad)', 'FontSize', 12);
ylabel('Lateral Force F_y (N)', 'FontSize', 12);
title(['Tire Curve Comparison (F_z = ' num2str(round(Fz_nominal)) ' N)'], 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 10);
xlim([-0.25, 0.25]); % Zoom in on the relevant range

% Add "Cornering Stiffness" annotation
% K = B*C*D*Fz
K_gt = prod(theta_gt(1:3)) * Fz_nominal;
K_fr = prod(theta_front(1:3)) * Fz_nominal;
K_re = prod(theta_rear(1:3)) * Fz_nominal;

dim = [0.15 0.6 0.3 0.3];
str = {['Stiffness (C_\alpha):'], ...
       ['GT:   ' num2str(round(K_gt))], ...
       ['Front: ' num2str(round(K_fr))], ...
       ['Rear:  ' num2str(round(K_re))]};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');

hold off;