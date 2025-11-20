function [x_smooth, P_smooth, P_cross_smooth] = RTS_smoother_EM(x_filt, P_filt, x_pred, P_pred, Ad_hist)
    [N, num_states] = size(x_filt);
    
    x_smooth = zeros(N, num_states);
    P_smooth = zeros(num_states, num_states, N);
    P_cross_smooth = zeros(num_states, num_states, N); % P_{t, t-1 | T}
    
    % Initialize with final state
    x_smooth(N, :) = x_filt(N, :);
    P_smooth(:, :, N) = P_filt(:, :, N);
    
    for k = (N - 1) : -1 : 1
        x_k_k = x_filt(k, :)';
        P_k_k = P_filt(:, :, k);
        
        x_k1_k = x_pred(k + 1, :)';
        P_k1_k = P_pred(:, :, k + 1);
        Ad_k = Ad_hist(:, :, k);

        jitter = 1e-9 * eye(size(P_k1_k));
        
        % Gain
        J = (P_k_k * Ad_k') / (P_k1_k + jitter);
        
        % State
        x_smooth(k, :) = (x_k_k + J * (x_smooth(k+1, :)' - x_k1_k))';
        
        % Covariance
        P_smooth(:, :, k) = P_k_k + J * (P_smooth(:, :, k+1) - P_k1_k) * J';
        
        % Lag-One Covariance (Approximation for EM)
        % P_{t+1, t | T} = P_{t+1 | T} * J_t'
        % We store this at index k+1 (representing relationship between k+1 and k)
        P_cross_smooth(:, :, k+1) = P_smooth(:, :, k+1) * J';
    end
end