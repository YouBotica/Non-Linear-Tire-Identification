% --- Plot Results Script ---
% Assumes 'track', 'sim_X', and 'sim_Y' are in your Workspace

figure;
hold on;

% 1. Plot the Track Path (Blue)
plot(track.X, track.Y, 'b-', 'LineWidth', 2, 'DisplayName', 'Track Path');

% 2. Plot the Vehicle's Path (Red, Dashed)
plot(out.Xpos, out.Ypos, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Vehicle Path');

% 3. Plot the Start Point
plot(out.Xpos(1), out.Ypos(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');

% 4. Make it look good
axis equal; % <-- Most important! Prevents a distorted oval.
grid on;
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Vehicle Path vs. Track');
legend('Location', 'best');
hold off;