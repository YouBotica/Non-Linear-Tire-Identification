function plot_EM_comparison(y_log, x_human_tuned, x_em_tuned, dt)
% Plots a comparison of States AND Disturbances before and after EM tuning.
%
% INPUTS:
%   y_log:          (N x 5) Matrix of raw measurements [vx, vy, r, FzF, FzR]
%   x_human_tuned:  (N x 6) Matrix of states from your initial guess
%   x_em_tuned:     (N x 6) Matrix of states from the final EM result
%   dt:             Sample time

    [N, ~] = size(y_log);
    time = (0:N-1) * dt;
    
    % Create figure with a 3x2 layout
    figure('Name', 'EM Tuning Comparison: States & Disturbances', 'Color', 'w', 'Position', [100 100 1200 800]);
    sgtitle('Impact of EM Tuning: Human (Red) vs. EM (Blue)', 'FontSize', 14);
    
    % --- ROW 1: LONGITUDINAL ---
    
    % Subplot 1: Vx State
    subplot(3, 2, 1);
    plot(time, y_log(:, 1), 'k.', 'MarkerSize', 4); hold on;
    plot(time, x_human_tuned(:, 1), 'r--', 'LineWidth', 1.5);
    plot(time, x_em_tuned(:, 1), 'b-', 'LineWidth', 1.5);
    ylabel('v_x (m/s)');
    title('Longitudinal Velocity');
    legend('Measured', 'Human-Tuned', 'EM-Tuned', 'Location', 'best');
    grid on;
    
    % Subplot 2: Xi_x Disturbance
    subplot(3, 2, 2);
    plot(time, x_human_tuned(:, 4), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, x_em_tuned(:, 4), 'b-', 'LineWidth', 1.5);
    ylabel('\xi_x (Disturbance)');
    title('Longitudinal Disturbance (\xi_x)');
    legend('Human-Tuned', 'EM-Tuned', 'Location', 'best');
    grid on;
    
    % --- ROW 2: LATERAL ---
    
    % Subplot 3: Vy State
    subplot(3, 2, 3);
    plot(time, y_log(:, 2), 'yo', 'MarkerSize', 4); hold on;
    plot(time, x_human_tuned(:, 2), 'r--', 'LineWidth', 1.5);
    plot(time, x_em_tuned(:, 2), 'b-', 'LineWidth', 1.5);
    ylabel('v_y (m/s)');
    title('Lateral Velocity');
    grid on;
    
    % Subplot 4: Xi_y Disturbance
    subplot(3, 2, 4);
    plot(time, x_human_tuned(:, 5), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, x_em_tuned(:, 5), 'b-', 'LineWidth', 1.5);
    ylabel('\xi_y (Disturbance)');
    title('Lateral Disturbance (\xi_y)');
    grid on;
    
    % --- ROW 3: YAW ---
    
    % Subplot 5: Yaw Rate State
    subplot(3, 2, 5);
    plot(time, y_log(:, 3), 'k.', 'MarkerSize', 4); hold on;
    plot(time, x_human_tuned(:, 3), 'r--', 'LineWidth', 1.5);
    plot(time, x_em_tuned(:, 3), 'b-', 'LineWidth', 1.5);
    ylabel('r (rad/s)');
    xlabel('Time (s)');
    title('Yaw Rate');
    grid on;
    
    % Subplot 6: Xi_psi Disturbance
    subplot(3, 2, 6);
    plot(time, x_human_tuned(:, 6), 'r--', 'LineWidth', 1.5); hold on;
    plot(time, x_em_tuned(:, 6), 'b-', 'LineWidth', 1.5);
    ylabel('\xi_\psi (Disturbance)');
    xlabel('Time (s)');
    title('Yaw Disturbance (\xi_\psi)');
    grid on;
end