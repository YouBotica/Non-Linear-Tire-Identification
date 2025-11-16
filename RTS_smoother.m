function [x_smooth, r_smooth] = RTS_smoother(logged_x_t_t, logged_r_t_t, ...
                                               logged_x_t_t_1, logged_r_t_t_1, ...
                                               logged_Ad)
% run_rts_smoother: Performs a full RTS backward pass.
%
% This function takes the *entire history* of data logged from a 
% Kalman Filter (your "Forward Pass") and computes the smoothed
% state and covariance estimates (the "Backward Pass").
%
% INPUTS:
%   logged_x_t_t:   (N x 6) matrix of updated states [x_t|t]
%   logged_r_t_t:   (6 x 6 x N) array of updated covariances [r_t|t]
%   logged_x_t_t_1: (N x 6) matrix of predicted states [x_t|t-1]
%   logged_r_t_t_1: (6 x 6 x N) array of predicted covariances [r_t|t-1]
%   logged_Ad:      (6 x 6 x N) array of discrete A matrices [Ad]
%
% OUTPUTS:
%   x_smooth:       (N x 6) matrix of smoothed states [x_t|T]
%   r_smooth:       (6 x 6 x N) array of smoothed covariances [r_t|T]

fprintf('Starting RTS Smoother Backward Pass...\n');

% Get dimensions
[N, num_states] = size(logged_x_t_t);

% --- 1. Pre-allocate Output Arrays ---
% We will store our smoothed estimates here.
x_smooth = zeros(N, num_states);
r_smooth = zeros(num_states, num_states, N);

% --- 2. Initialize the Backward Pass ---
% The smoothed state at the final time step (T) is just
% the final updated state from the Kalman Filter.
x_smooth(N, :) = logged_x_t_t(N, :);
r_smooth(:, :, N) = logged_r_t_t(:, :, N);

% --- 3. Run the Backward Loop ---
% Iterate from k = T-1 down to 1.
% (Note: k in our loop corresponds to t-1 in your book's notation)
for k = (N - 1) : -1 : 1
    
    % Get the logged data for time k (which is t-1)
    x_k_k = logged_x_t_t(k, :)'; % (6x1) column vector
    r_k_k = logged_r_t_t(:, :, k); % (6x6) matrix
    
    % Get the logged data for time k+1 (which is t)
    x_k1_k = logged_x_t_t_1(k + 1, :)'; % (6x1) column vector
    r_k1_k = logged_r_t_t_1(:, :, k + 1); % (6x6) matrix
    Ad_k   = logged_Ad(:, :, k);          % (6x6) matrix
    
    % Get the *previously smoothed* data from k+1
    x_smooth_k1 = x_smooth(k + 1, :)'; % (6x1) column vector
    r_smooth_k1 = r_smooth(:, :, k + 1); % (6x6) matrix

    % --- 4. Calculate Smoother Gain (J) ---
    % (This is your "Future Conditioning" step)
    % We use the more numerically stable / (right-divide) instead of inv()
    % J_k = r_k_k * Ad_k' * inv(r_k1_k);
    J_k = (r_k_k * Ad_k') / r_k1_k;

    % --- 5. Calculate Smoothed State & Covariance ---
    % (This is your "Backward Step")
    
    % Smoothed state: x_{t-1|T} = x_{t-1|t-1} + J_{t-1}(x_{t|T} - x_{t|t-1})
    x_smooth_k = x_k_k + J_k * (x_smooth_k1 - x_k1_k);
    
    % Smoothed covariance: r_{t-1|T} = r_{t-1|t-1} + J_{t-1}(r_{t|T} - r_{t|t-1})J_{t-1}'
    r_smooth_k = r_k_k + J_k * (r_smooth_k1 - r_k1_k) * J_k';
    
    % Store the results (transpose state back to a row)
    x_smooth(k, :) = x_smooth_k';
    r_smooth(:, :, k) = r_smooth_k;
    
end

fprintf('RTS Smoother finished.\n');
end