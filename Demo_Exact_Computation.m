% Written by: Mehryar Emambakhsh
% Email: mehryar_emam@yahoo.com
% Date: 27 Sep 2017
% Paper:
% P. Garcia, M. Emambakhsh, A. Wallace, “Learning to Approximate Computing at Run-time,”
% in IET 3rd International Conference on Intelligent Signal Processing (ISP
% 2017), 2017, to appear.
function Demo_Exact_Computation
% This code simulates circular rotation for a target which is being tracked
% by an extended Kalman filter. No approximation is performed here, i.e
% exact computation.
warning off
clc
% close all

% It is assumed that the sensor is located at the origin with FoV of XY
% plane (x > 0 and y > 0).

%% Establish the scene and initialisation
% Target's circular motion parameters
my_radius = 6;
centre_of_rotation = [7.5, 7.5];
initial_target_location = [10, 5];
delta_theta = 1* pi/180;

% Initialisations for target
theta_true_t_minus_1 = 0;
x_true_t_minus_1 = initial_target_location(1);
y_true_t_minus_1 = initial_target_location(2);

% ranges initialisations
min_x = 0;
max_x = 15;
min_y = 0;
max_y = 15;

measurement_noise_intensity = 0.1;

%%%%%%%%%%%%%%%%%%% Plot
figure; %plot(initial_target_location(1), initial_target_location(2), 'r.', 'markersize', 25)
xlim([min_x, max_x]), ylim([min_y, max_y])
xlabel('x')
ylabel('y')
% title('Red: estimated state, Blue: true state')
%%%%%%%%%%%%%%%%%%%%%%%%


%% Main loop

Concat_to_plot = [];
% while(1)
for iter_cnt = 1: 365
    %%%%%%%%%%%%%%%%%%% Update the scene, with the new target's lcoation
    x_true_t = my_radius* cos(theta_true_t_minus_1) + centre_of_rotation(1);
    y_true_t = my_radius* sin(theta_true_t_minus_1) + centre_of_rotation(2);
    theta_true_t = theta_true_t_minus_1 + delta_theta;
    X_true_t = [x_true_t; y_true_t];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%% Sensor measurement: Obtain the measurement: Z_t
    % Note that for the simulation, we are using the true state of the
    % target's motion, with an added Gaussian noise:
    r_k = sqrt(x_true_t^ 2 + y_true_t^ 2);
    phi_k = atan2(y_true_t, x_true_t);
    measurement_variance = measurement_noise_intensity* rand(2, 1);
    m_k = randn(2, 1) .* measurement_variance;
    %     m_k = 0;
    Z_k = [r_k; phi_k] + m_k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%% Extended Kalman filter step
    % Inputs are: Z_k
    [X_k_exact, ~] = myEKFEstimator(Z_k);
    %%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% Update target's angular movement
    theta_true_t_minus_1 = theta_true_t;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% Plotting step
    if mod(iter_cnt, 1) == 0
        Concat_to_plot = [Concat_to_plot, X_true_t];
        EndToPlot = 40;
        if size(Concat_to_plot, 2) < EndToPlot
            plot(Concat_to_plot(1, :), Concat_to_plot(2, :), '.r', 'markersize', 18),
        else
            plot(Concat_to_plot(1, end-EndToPlot+1:end), Concat_to_plot(2, end-EndToPlot+1:end), '.r', 'markersize', 18),
        end
        hold on
        plot(X_k_exact(1), X_k_exact(2), '.g', 'markersize', 18),
        
        xlim([min_x, max_x]), ylim([min_y, max_y])
        drawnow
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
title('Red: true state tragectory; Green: EKF tracking')
end
