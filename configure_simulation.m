% --- 1. Clear workspace and close figures ---
clear; clc; close all;

% --- 2. Define Track Parameters ---
longStraight = 1200;   % meters
shortStraight = 500;  % meters
turnRadius = 250;      % meters
pointsPerSegment = 500; % 100 points for each of the 8 segments
v_straight = 20;
v_turn = 10;

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