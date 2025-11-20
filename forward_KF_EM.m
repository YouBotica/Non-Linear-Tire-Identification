function [x_hist, P_hist, x_pred_hist, P_pred_hist, Ad_hist] = forward_KF_EM(y, u, Q, R, c, dt)
    [N, ~] = size(y);
    num_states = 6;
    
    % Initialize
    x = zeros(num_states, 1);
    x(1) = 5.0; % Initial velocity guess
    P = eye(num_states) * 1.0;
    
    % Storage
    x_hist = zeros(N, num_states);
    P_hist = zeros(num_states, num_states, N);
    x_pred_hist = zeros(N, num_states);
    P_pred_hist = zeros(num_states, num_states, N);
    Ad_hist = zeros(num_states, num_states, N);
    
    C = [eye(3), zeros(3,3)];
    
    for k = 1:N
        % Get inputs for this step
        uk = u(k, :)';
        yk = y(k, 1:3)'; % Only measure vx, vy, r
        Fz = y(k, 4:5);  % FzF, FzR (parameters, not states)
        
        % --- Linearize System (Calculate Ad, Bd) ---
        [Ad, Bd] = get_linear_model_internal(x, uk, Fz, c, dt);
        
        % --- Predict ---
        x_pred = Ad * x + Bd * uk;
        P_pred = Ad * P * Ad' + Q;
        
        % --- Update ---
        % Joseph form for stability
        K = P_pred * C' / (C * P_pred * C' + R);
        x = x_pred + K * (yk - C * x_pred);
        I = eye(num_states);
        P = (I - K * C) * P_pred * (I - K * C)' + K * R * K';
        
        % --- Store ---
        x_hist(k, :) = x';
        P_hist(:, :, k) = P;
        x_pred_hist(k, :) = x_pred';
        P_pred_hist(:, :, k) = P_pred;
        Ad_hist(:, :, k) = Ad;
    end
end

function [Ad, Bd] = get_linear_model_internal(x, u, Fz, c, dt)
    % Same linear model logic as in main script
    vx = max(x(1), 1.0);
    FzF = Fz(1); FzR = Fz(2);
    Caf = c.K_pacejka * FzF; Car = c.K_pacejka * FzR;
    
    Yv = -(2*Caf + 2*Car) / (c.m * vx);
    Yr = -vx - (2*Caf*c.lf - 2*Car*c.lr) / (c.m * vx);
    Nv = -(2*Caf*c.lf - 2*Car*c.lr) / (c.Izz * vx);
    Nr = -(2*Caf*c.lf^2 + 2*Car*c.lr^2) / (c.Izz * vx);
    
    A = [ -c.Cd/c.m, 0, 0, 1/c.m, 0, 0;
          0, Yv, Yr, 0, 1/c.m, 0;
          0, Nv, Nr, 0, 0, 1/c.Izz;
          zeros(3, 6) ];
    B = [ 1/c.m, 0; 0, 2*Caf/c.m; 0, 2*Caf*c.lf/c.Izz; zeros(3, 2) ];
    
    Ad = eye(6) + A*dt;
    Bd = B*dt;
end