% Written by: Mehryar Emambakhsh
% Email: mehryar_emam@yahoo.com
% Date: 27 Sep 2017
% Paper:
% P. Garcia, M. Emambakhsh, A. Wallace, “Learning to Approximate Computing at Run-time,”
% in IET 3rd International Conference on Intelligent Signal Processing (ISP
% 2017), 2017, to appear.
function [all_ERR, all_KL, all_energy, app_level_hist] = Dynamic_approximation_detailed_visualisation(my_KL_thresh)
% This function performs dynamic approximation using the given KL threshold
% as input.
% Compared to Dynamic_approximation.m, it provides more detailed 
% visualisation of the compromisation between the approximation level, 
% power consumtion and tracking error.
% my_KL_thresh: 1 X 1 scalar
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
figure; subplot(1, 2, 1) %plot(initial_target_location(1), initial_target_location(2), 'r.', 'markersize', 25)
xlim([min_x, max_x]), ylim([min_y, max_y])
xlabel('x')
ylabel('y')
% title('Red: estimated state, Blue: true state')
%%%%%%%%%%%%%%%%%%%%%%%%

% Number of stacked states to compute histogram for KL Div
N_stack = 75;
mu_ref = my_radius;
var_ref = 0.25;

% Defining levels of approximation : the higher index in the cell, the
% more intense approximation
app_levels = {'E', 'D', 'C', 'ACDE' , 'ACD' , 'AD' , 'AC' , 'A',...
    'AE', 'ADE', 'ACE', 'ACF', 'ACEF', 'AEF', 'ACDF', 'ACDEF', ...
    'AF', 'ADEF', 'ADF'};
app_levels = app_levels(1: end-6);

% Computed energy from FPGA implementation
energy_vector = [0.781, 0.77, 0.77, .887359198999* 0.799,...
.927409261577* 0.799,.946182728411* 0.799,.946182728411* 0.799,.982478097622* 0.799,...
.959949937422* 0.799,.92365456821* 0.799,.959949937422* 0.799,.921151439299* 0.799,...
.898623279099* 0.799,.93491864831* 0.799,.884856070088* 0.799,.862327909887* 0.799,...
.957446808511* 0.799,.898623279099* 0.799,.921151439299* 0.799];

% my_KL_thresh = 4;
curr_app_level = 1;

%% Main loop

Concat_to_plot = [];
all_ERR = [];
all_KL = [];
stacked_samples = [];
all_energy = [];
app_level_hist = zeros(1, 19);
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
    % The previous state vector is defined as persistent within the
    % function (persistent is static in Matlab!)
    [X_k, cov_X_new] = Approximate_Kalman(Z_k, app_levels{curr_app_level});
    [X_k_exact, ~] = myEKFEstimator(Z_k);
    %%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% Update target's angular movement
    theta_true_t_minus_1 = theta_true_t;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% Plotting step
    
    if mod(iter_cnt, 1) == 0
        subplot(1, 2, 1)
        Concat_to_plot = [Concat_to_plot, X_true_t];
        EndToPlot = 40;
        if size(Concat_to_plot, 2) < EndToPlot
            plot(Concat_to_plot(1, :), Concat_to_plot(2, :), '.r', 'markersize', 18),
        else
            plot(Concat_to_plot(1, end-EndToPlot+1:end), Concat_to_plot(2, end-EndToPlot+1:end), '.r', 'markersize', 18),
        end
        hold on
%         plot(Concat_to_plot(1, end), Concat_to_plot(2, end), '*g', 'markersize', 22)
        plot(X_k(1), X_k(2), '.b', 'markersize', 18), 
        plot(X_k_exact(1), X_k_exact(2), '.g', 'markersize', 18), 
        
%         hold off
        % Plotting ellipses
%         ell_points = plott_ellipse(X_k, cov_X_new);
%         hold on
%         plot(ell_points(1, :), ell_points(2, :), 'r')
%         hold off
        xlim([min_x, max_x]), ylim([min_y, max_y])
        drawnow
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%% Computing and concatanating  accuracy error in real-time
    current_radius_of_states = sqrt((X_k(1) - centre_of_rotation(1))^2 + (X_k(2) - centre_of_rotation(2))^2);
    stacked_samples = [stacked_samples, current_radius_of_states];
    if length(stacked_samples) >= N_stack
        curr_stacked_samples = stacked_samples(end: -1: end - 50 + 1);
        curr_window_mean = mean(curr_stacked_samples);
        curr_window_var = var(curr_stacked_samples);
        curr_KL_output = myKLDiv_computation(mu_ref, var_ref, curr_window_mean, curr_window_var);
        title(['Current KL value: ' num2str(curr_KL_output)]), drawnow
        
        %%%%%%%%%%%%%% Compute the consumed energy: assuming delta_t = 0.1;
        all_energy = [all_energy, energy_vector(curr_app_level)* 0.1];
        app_level_hist(curr_app_level) = app_level_hist(curr_app_level) + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%% Now approximate based on the value of KL div:
        if curr_KL_output > my_KL_thresh
            % If KL-div was high, it means the distributions are so dissimilar,
            % so do less intense approximation
            if curr_app_level ~= 1
                curr_app_level = curr_app_level - 1;
            else
                curr_app_level = 1;
            end
        else
            % We are cool! Do more approximation
            if curr_app_level ~= length(app_levels)
                curr_app_level = curr_app_level + 1;
            else
                curr_app_level = length(app_levels);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        myErr = sqrt(sum(([x_true_t; y_true_t] - X_k(1: 2)).^2));
        all_ERR = [all_ERR, myErr];
        all_KL = [all_KL, curr_KL_output];
        
        %%%%% Plotting the energy, error and approximation level as bar
        %%%%% graph
        subplot(1, 2, 2)
        bar1 = bar(1, 100* curr_app_level/ length(app_levels));
        hold on
        bar2 = bar(2, 100* energy_vector(curr_app_level)/max(energy_vector));
        hold off
        set(bar1, 'FaceColor', 'green');
        set(bar2, 'FaceColor', 'red');
        
        set(gca, 'XTick', [1, 2])
        set(gca, 'XTickLabel', {'Approximation level', 'Power'})
        
        ylabel('Percent (%)')
        ylim([1, 100])
        legend('Approximation level', 'Power')
        title(sprintf('Normalised Approximation level and \n its corrsponding power consumtion'))
        drawnow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
title('Red: true state tragectory; Green: EKF tracking; Blue: Approximate computing result')
end

function curr_KL_output = myKLDiv_computation(mu_ref, var_ref, mu_current, var_current)
curr_KL_output = log10(var_ref/(eps+double(var_current))) + (var_current^2 + (mu_current - mu_ref)^2)/(2* var_ref^2) - 0.5;
end


function ell_points = plott_ellipse(X_k, cov_X_new)
xCenter = double(X_k(1)); yCenter = double(X_k(2));
[D, V] = eig(double(cov_X_new(1: 2, 1: 2)));
major_axis = D(:, 1);
xRadius = sqrt(V(1));
yRadius = sqrt(V(4));
mytheta = 0 : 0.01 : 2*pi;
x = xRadius * cos(mytheta) + xCenter;
y = yRadius * sin(mytheta) + yCenter;

my_tilt_angle = atan2(major_axis(2), major_axis(1));

mean_matrix = repmat(mean([x; y], 2), 1, length(x));
no_mean = [x; y] - mean_matrix;
ell_points = [cos(my_tilt_angle) -sin(my_tilt_angle); sin(my_tilt_angle) cos(my_tilt_angle)]* no_mean;
ell_points = ell_points + mean_matrix;

end