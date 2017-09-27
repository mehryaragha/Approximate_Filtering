% Written by: Mehryar Emambakhsh
% Email: mehryar_emam@yahoo.com
% Date: 27 Sep 2017
% Paper:
% P. Garcia, M. Emambakhsh, A. Wallace, “Learning to Approximate Computing at Run-time,”
% in IET 3rd International Conference on Intelligent Signal Processing (ISP
% 2017), 2017, to appear.
function [X_k, cov_X_new] = myEKFEstimator(Z_k)
% This function performs Extended Kalman filter: prediction is performed
% over the current state (t-1), which is defined as static (persistent)
% variable within the function. The update is performed using the current
% mesuarement Z_k.

persistent X_k_t_minus_1 cov_X_k_minus_1
% At the first iteration the static variable is empty, therefore use Z_k to
% find its value:
if isempty(X_k_t_minus_1)
    r_k = Z_k(1);
    phi_k = Z_k(2);
    X_k_t_minus_1 = [r_k* cos(phi_k); r_k* sin(phi_k); pi/2];
    cov_X_k_minus_1 = eye(3);
end

x_k_minus_1 = X_k_t_minus_1(1);
y_k_minus_1 = X_k_t_minus_1(2);
theta_k_minus_1 = X_k_t_minus_1(3);

v_k = 1; % radial speed in m/sec
w_k = 10* pi/180; % angular speed in rad/sec

%%%%%%%%%%% Some intialisation
% Estimation's covariance matrix
delta_t = 0.1; % time resolution
% Q_k_minus_1 = eye(3)* state_noise_intensity* rand* abs(randn);
Q_k_minus_1 = 2*eye(3);
% Measurement noise covariance matrix
% R_k = eye(2)* measurement_noise_intensity* rand* abs(randn);
R_k = 2*eye(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Predict the next state
state_noise_intensity = 0;
state_noise_variance = state_noise_intensity* rand(3, 1);
% Because the third element of state_noise_variance corresponds to
% the robot pose, I assume the given intensity was in degree. So here
% before I proceed, it is changed to radian.
state_noise_variance(3) = state_noise_variance(3)* pi/ 180;
n_k = randn(3, 1) .* state_noise_variance;

% n_k = 0;
x_k = x_k_minus_1 + delta_t* v_k* cos(w_k* delta_t + theta_k_minus_1);
y_k = y_k_minus_1 + delta_t* v_k* sin(w_k* delta_t + theta_k_minus_1);
theta_k = theta_k_minus_1 + w_k* delta_t;
X_k = [x_k; y_k; theta_k] + n_k;
% Compute the estimation's covariance via Jacobian computation.
J_X_k_minus_1 = eye(3) +...
    [0 0 -delta_t* v_k* sin(w_k* delta_t + theta_k_minus_1);
    0 0 delta_t* v_k* cos(w_k* delta_t + theta_k_minus_1)
    0 0 0];
cov_X_k = J_X_k_minus_1* cov_X_k_minus_1* J_X_k_minus_1' + Q_k_minus_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Correction step
% Compute Jacobian
H_k = (1/ sqrt(x_k^2 + y_k^2)) * ...
    [x_k, y_k, 0;
    -y_k, x_k, 0];

% Compute Kalman gain
S_k = H_k* cov_X_k* H_k' + R_k;
K_k = cov_X_k* H_k'/ (S_k);

X_k_new = X_k + K_k* (Z_k - [sqrt(x_k^ 2 + y_k^ 2); atan2(y_k, x_k)]);
cov_X_new = (eye(3) - K_k* H_k)* cov_X_k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Update the state
X_k_t_minus_1 = X_k_new;
cov_X_k_minus_1 = cov_X_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end