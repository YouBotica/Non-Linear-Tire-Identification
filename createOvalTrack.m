% Save this as 'createRovalTrack_with_speed.m'
function track = createOvalTrack(longStraight, shortStraight, turnRadius, ...
                                           pointsPerSegment, v_straight, v_turn)
    % Creates a "roval" track with a pre-calculated, smooth speed profile.
    % This version does NOT use friction 'mu'. It uses the speeds you give it.
    %
    % Inputs:
    %   longStraight, shortStraight, turnRadius, pointsPerSegment: Track geometry
    %   v_straight: The target speed for the long straights (e.g., 20 m/s)
    %   v_turn:     The target speed for the turns (e.g., 8 m/s)
    
    % --- 1. Generate Track Geometry (Same as before) ---
    y_top = shortStraight/2 + turnRadius;
    x1 = linspace(0, longStraight, pointsPerSegment)';
    y1 = ones(pointsPerSegment, 1) * y_top;
    yaw1 = zeros(pointsPerSegment, 1);
    t = linspace(pi/2, 0, pointsPerSegment)';
    center_x = longStraight; center_y = shortStraight/2;
    x2 = center_x + turnRadius * cos(t);
    y2 = center_y + turnRadius * sin(t);
    yaw2 = -((pi/2) - t);
    x_right = longStraight + turnRadius;
    x3 = ones(pointsPerSegment, 1) * x_right;
    y3 = linspace(shortStraight/2, -shortStraight/2, pointsPerSegment)';
    yaw3 = ones(pointsPerSegment, 1) * (-pi/2);
    t = linspace(0, -pi/2, pointsPerSegment)';
    center_x = longStraight; center_y = -shortStraight/2;
    x4 = center_x + turnRadius * cos(t);
    y4 = center_y + turnRadius * sin(t);
    yaw4 = t - pi/2;
    y_bottom = -shortStraight/2 - turnRadius;
    x5 = linspace(longStraight, 0, pointsPerSegment)';
    y5 = ones(pointsPerSegment, 1) * y_bottom;
    yaw5 = ones(pointsPerSegment, 1) * pi;
    t = linspace(-pi/2, -pi, pointsPerSegment)';
    center_x = 0; center_y = -shortStraight/2;
    x6 = center_x + turnRadius * cos(t);
    y6 = center_y + turnRadius * sin(t);
    yaw6 = t - pi/2;
    x_left = -turnRadius;
    x7 = ones(pointsPerSegment, 1) * x_left;
    y7 = linspace(-shortStraight/2, shortStraight/2, pointsPerSegment)';
    yaw7 = ones(pointsPerSegment, 1) * (pi/2);
    t = linspace(pi, pi/2, pointsPerSegment)';
    center_x = 0; center_y = shortStraight/2;
    x8 = center_x + turnRadius * cos(t);
    y8 = center_y + turnRadius * sin(t);
    yaw8 = (t - pi) + pi/2;
    
    track.X = [x1; x2; x3; x4; x5; x6; x7; x8];
    track.Y = [y1; y2; y3; y4; y5; y6; y7; y8];
    track.Yaw = [yaw1; yaw2; yaw3; yaw4; yaw5; yaw6; yaw7; yaw8];
    track.numPoints = length(track.X);
    
    % --- 2. Generate Speed Profile ---
    fprintf('Generating speed profile...\n');
    fprintf('  Straight speed: %.1f m/s\n', v_straight);
    fprintf('  Turn speed: %.1f m/s\n', v_turn);
    
    % Define acceleration/braking (e.g., 2.0 m/s^2)
    max_accel = 5.0; 
    
    % Create a "base" speed profile
    % 1=straight, 2=turn, 3=straight, 4=turn, etc.
    segment_speeds = [v_straight; v_turn; v_straight; v_turn; ...
                      v_straight; v_turn; v_straight; v_turn];
                  
    speed_profile = zeros(track.numPoints, 1);
    for i = 1:8
        start_idx = (i-1)*pointsPerSegment + 1;
        end_idx = i*pointsPerSegment;
        speed_profile(start_idx:end_idx) = segment_speeds(i);
    end
    
    % --- 3. Smooth the profile with braking/acceleration ---
    % (This logic is good and remains the same)
    
    % First pass (braking): iterate backwards
    for i = (track.numPoints - 1):-1:1
        ds = sqrt((track.X(i+1)-track.X(i))^2 + (track.Y(i+1)-track.Y(i))^2);
        v_max_braking = sqrt(speed_profile(i+1)^2 + 2 * max_accel * ds);
        speed_profile(i) = min(speed_profile(i), v_max_braking);
    end
    
    % Second pass (acceleration): iterate forwards
    for i = 1:(track.numPoints - 1)
        ds = sqrt((track.X(i+1)-track.X(i))^2 + (track.Y(i+1)-track.Y(i))^2);
        v_max_accel = sqrt(speed_profile(i)^2 + 2 * max_accel * ds);
        speed_profile(i+1) = min(speed_profile(i+1), v_max_accel);
    end
    
    speed_profile(end) = speed_profile(1); % Handle wrap-around
    
    track.SpeedProfile = speed_profile;
    fprintf('Speed profile generated.\n');
end