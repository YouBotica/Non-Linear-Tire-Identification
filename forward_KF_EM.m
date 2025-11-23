function [x_hist, P_hist, x_pred_hist, P_pred_hist, Ad_hist] = forward_KF_EM(y, u, Q, R, P_initial, c, dt)
    [N, ~] = size(y);
    num_states = 6;
    
    % Initialize
    x = zeros(num_states, 1);
    x(1) = 5.0; % Initial velocity guess
    P = P_initial;
    
    % Storage
    x_hist = zeros(N, num_states);
    P_hist = zeros(num_states, num_states, N);
    x_pred_hist = zeros(N, num_states);
    P_pred_hist = zeros(num_states, num_states, N);
    Ad_hist = zeros(num_states, num_states, N);
    
    for k = 1:N
        % Get inputs for this step
        uk = u(k, :)';
        yk = y(k, 1:4)'; % Only measure vx, vy, r, ay
        Fz = y(k, 5:6);  % FzF, FzR (parameters, not states)
        
        % --- Linearize System (Calculate Ad, Bd) ---
        [Ad, Bd, Yv, Yr] = get_linear_model_internal(x, uk, Fz, c, dt);
        

        C_aug = [ % C matrix changes over time, more specifically its last row.
           eye(3), zeros(3);
           0, Yv, (Yr + x(1)),  0, 1/c.m, 0;    % Computes ay (NEW)
        ];

        % --- Predict ---
        x_pred = Ad * x + Bd * uk;
        P_pred = Ad * P * Ad' + Q;
        
        % --- Update ---
        % Joseph form for stability
        K = P_pred * C_aug' / (C_aug * P_pred * C_aug' + R);
        x = x_pred + K * (yk - C_aug * x_pred);
        I = eye(num_states);
        P = (I - K * C_aug) * P_pred * (I - K * C_aug)' + K * R * K';
        
        % --- Store ---
        x_hist(k, :) = x';
        P_hist(:, :, k) = P;
        x_pred_hist(k, :) = x_pred';
        P_pred_hist(:, :, k) = P_pred;
        Ad_hist(:, :, k) = Ad;
    end
end

function [Ad, Bd, Yv, Yr] = get_linear_model_internal(x, u, Fz, c, dt)
    % Same linear model logic as in main script
    vx = max(x(1), 1.0);
    FzF = Fz(1); FzR = Fz(2);
    Caf = c.K_pacejka * FzF; Car = c.K_pacejka * FzR;
    
    Yv = -(Caf + Car) / (c.m * vx);
    Yr = -vx - (Caf*c.lf - Car*c.lr) / (c.m * vx);
    Nv = -(Caf*c.lf - Car*c.lr) / (c.Izz * vx);
    Nr = -(Caf*c.lf^2 + Car*c.lr^2) / (c.Izz * vx);
    
    A = [ -c.Cd/c.m, 0, 0, 1/c.m, 0, 0;
          0, Yv, Yr, 0, 1/c.m, 0;
          0, Nv, Nr, 0, 0, 1/c.Izz;
          zeros(3, 6) ];
    B = [ 1/c.m, 0; 0, Caf/c.m; 0, Caf*c.lf/c.Izz; zeros(3, 2) ];
    
    % Discretize:
    Ad = eye(6) + A*dt;
    Bd = B*dt;
end