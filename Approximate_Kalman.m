% Written by: Mehryar Emambakhsh
% Email: mehryar_emam@yahoo.com
% Date: 27 Sep 2017
% Paper:
% P. Garcia, M. Emambakhsh, A. Wallace, “Learning to Approximate Computing at Run-time,”
% in IET 3rd International Conference on Intelligent Signal Processing (ISP
% 2017), 2017, to appear.
function [X_k, cov_X_new] = Approximate_Kalman(Z_k, app_level)
% This function is now capable of using Matlab's fi function to change the
% floating point precision to perform approximate computing.
% Also, the app_level is now given as a strings: (Namings according to Paulo's)
% app_level = 'A' means do the fi Bit Width 27 15
% app_level = 'B' means do the fi Bit Width 18 1
% app_level = 'C' approximate cosine
% app_level = 'D' approximate sine
% app_level = 'E' approximate SQRT
% app_level = 'F' approximate ATAN
% Combination is allowed, e.g. if app_level = 'AD', means do
% precision 'A' and approximate sine only.

if nargin == 1
    app_level = '';
end

if ~isempty(strfind(app_level, 'C'))
    % Approximate cosine
    doCOSINE = true;
else
    doCOSINE = false;
end

if ~isempty(strfind(app_level, 'D'))
    % Approximate cosine
    doSINE = true;
else
    doSINE = false;
end

if ~isempty(strfind(app_level, 'E'))
    % Approximate cosine
    doSQRT = true;
else
    doSQRT = false;
end

if ~isempty(strfind(app_level, 'F'))
    % Approximate cosine
    doATAN = true;
else
    doATAN = false;
end

persistent X_k_t_minus_1 cov_X_k_minus_1

if ~isempty(strfind(app_level, 'B'))
    % Do precision 'A' explained above
    X_k_t_minus_1 = fi(X_k_t_minus_1, 1, 18, 17);
    cov_X_k_minus_1 = fi(cov_X_k_minus_1, 1, 18, 17);
    Z_k = fi(Z_k, 1, 18, 17);
end

if ~isempty(strfind(app_level, 'A'))
    % Do precision 'B' explained above
    X_k_t_minus_1 = fi(X_k_t_minus_1, 1, 27, 12);
    cov_X_k_minus_1 = fi(cov_X_k_minus_1, 1, 27, 12);
    Z_k = fi(Z_k, 1, 27, 12);
end


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
state_noise_intensity = 0.01;

% Here we assume the noise vector for estimation is fixed:
state_noise_variance = state_noise_intensity* ones(3, 1);
% Because the third element of state_noise_variance corresponds to
% the robot pose, I assume the given intensity was in degree. So here
% before I proceed, it is changed to radian.
state_noise_variance(3) = state_noise_variance(3)* pi/ 180;

n_k = ones(3, 1) .* state_noise_variance;

%APPROXIMATION 3
wrapped_input_angle = mywrapTo2Pi(w_k* delta_t + theta_k_minus_1, app_level);
if(doCOSINE)
    %let's replace the COS with a look up table
    myCos = getCos_FromLookUp(wrapped_input_angle);
else
    myCos = cos(w_k* delta_t + theta_k_minus_1);
end
x_k = x_k_minus_1 + delta_t* v_k* myCos;

%APPROXIMATION 4
if(doSINE)
    %let's replace the SIN with a look up table
    mySin = getSin_FromLookUp(wrapped_input_angle);
else
    mySin = sin(w_k* delta_t + theta_k_minus_1);
end
y_k = y_k_minus_1 + delta_t* v_k* mySin;


theta_k = theta_k_minus_1 + w_k* delta_t;
X_k = [x_k; y_k; theta_k] + n_k;
% Compute the estimation's covariance via Jacobian computation.


%APPROXIMATION 5
J_X_k_minus_1 = eye(3) +...
    [0 0 -delta_t* v_k* mySin;
    0 0 delta_t* v_k* myCos;
    0 0 0];

cov_X_k = J_X_k_minus_1* cov_X_k_minus_1* J_X_k_minus_1' + Q_k_minus_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Correction step
% Compute Jacobian



%APPROXIMATION 6
if(doSQRT)
    %let's replace the SQRT with a sum of absolutes
    H_k = (1/ (abs(x_k) + abs(y_k))) * ...
        [x_k, y_k, 0;
        -y_k, x_k, 0];
else
    H_k = (1/ sqrt(x_k^2 + y_k^2)) * ...
        [x_k, y_k, 0;
        -y_k, x_k, 0];
end

% Compute Kalman gain
S_k = H_k* cov_X_k* H_k' + R_k;
K_k = double(cov_X_k)* double(H_k)'/ double(S_k);

if ~isempty(strfind(app_level, 'B'))
    % Do precision 'A' explained above
    K_k = fi(K_k, 1, 18, 17);
end

if ~isempty(strfind(app_level, 'A'))
    % Do precision 'B' explained above
    K_k = fi(K_k, 1, 27, 12);
end


%APPROXIMATION 7
if(doATAN)
    %let's replace the SQRT with a sum of absolutes and ATAN2 with a look up table
    % Mehryar's comment: for atan, it's easier to approximate via taylor
    % series:
    
    % Checking if myAtan2 approximated here matches matlab's atan2, which
    % is computed from the positive side of x-axis:
    y_k_abs = abs(y_k);
    x_k_abs = abs(x_k);
    myAtan2 = y_k_abs/(eps+x_k_abs) - (1/3)*(y_k_abs/(eps+x_k_abs))^3 + (1/5)*(y_k_abs/(eps+x_k_abs))^5;
    if (y_k > 0 && x_k < 0)
        myAtan2 = pi - myAtan2;
    elseif (y_k < 0 && x_k < 0)
        myAtan2 = 2*pi + myAtan2;
    end
else
    myAtan2 = atan2(y_k, x_k);
end
X_k_new = X_k + K_k* (Z_k - [sqrt(x_k^ 2 + y_k^ 2); myAtan2]);

cov_X_new = (eye(3) - K_k* H_k)* cov_X_k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Update the state
X_k_t_minus_1 = double(X_k_new);
cov_X_k_minus_1 = double(cov_X_new);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function myCos = getCos_FromLookUp(wrapped_input_angle)
CosTable = [1.0000    0.9994    0.9976    0.9945    0.9903    0.9848    0.9781    0.9703    0.9613    0.9511    0.9397    0.9272    0.9135    0.8988...
    0.8829    0.8660    0.8480    0.8290    0.8090    0.7880    0.7660    0.7431    0.7193    0.6947    0.6691    0.6428    0.6157    0.5878...
    0.5592    0.5299    0.5000    0.4695    0.4384    0.4067    0.3746    0.3420    0.3090    0.2756    0.2419    0.2079    0.1736    0.1392...
    0.1045    0.0698    0.0349    0.0000   -0.0349   -0.0698   -0.1045   -0.1392   -0.1736   -0.2079   -0.2419   -0.2756   -0.3090   -0.3420...
    -0.3746   -0.4067   -0.4384   -0.4695   -0.5000   -0.5299   -0.5592   -0.5878   -0.6157   -0.6428   -0.6691   -0.6947   -0.7193   -0.7431...
    -0.7660   -0.7880   -0.8090   -0.8290   -0.8480   -0.8660   -0.8829   -0.8988   -0.9135   -0.9272   -0.9397   -0.9511   -0.9613   -0.9703...
    -0.9781   -0.9848   -0.9903   -0.9945   -0.9976   -0.9994   -1.0000   -0.9994   -0.9976   -0.9945   -0.9903   -0.9848   -0.9781   -0.9703...
    -0.9613   -0.9511   -0.9397   -0.9272   -0.9135   -0.8988   -0.8829   -0.8660   -0.8480   -0.8290   -0.8090   -0.7880   -0.7660   -0.7431...
    -0.7193   -0.6947   -0.6691   -0.6428   -0.6157   -0.5878   -0.5592   -0.5299   -0.5000   -0.4695   -0.4384   -0.4067   -0.3746   -0.3420...
    -0.3090   -0.2756   -0.2419   -0.2079   -0.1736   -0.1392   -0.1045   -0.0698   -0.0349   -0.0000    0.0349    0.0698    0.1045    0.1392...
    0.1736    0.2079    0.2419    0.2756    0.3090    0.3420    0.3746    0.4067    0.4384    0.4695    0.5000    0.5299    0.5592    0.5878...
    0.6157    0.6428    0.6691    0.6947    0.7193    0.7431    0.7660    0.7880    0.8090    0.8290    0.8480    0.8660    0.8829    0.8988...
    0.9135    0.9272    0.9397    0.9511    0.9613    0.9703    0.9781    0.9848    0.9903    0.9945    0.9976    0.9994    1.0000];

myInd = theta_to_ind_for_Cos_and_Sin(wrapped_input_angle, length(CosTable));
myCos = CosTable(myInd);
end

function mySin = getSin_FromLookUp(wrapped_input_angle)
SinTable = [0    0.0349    0.0698    0.1045    0.1392    0.1736    0.2079    0.2419    0.2756    0.3090    0.3420    0.3746    0.4067    0.4384...
    0.4695    0.5000    0.5299    0.5592    0.5878    0.6157    0.6428    0.6691    0.6947    0.7193    0.7431    0.7660    0.7880    0.8090...
    0.8290    0.8480    0.8660    0.8829    0.8988    0.9135    0.9272    0.9397    0.9511    0.9613    0.9703    0.9781    0.9848    0.9903...
    0.9945    0.9976    0.9994    1.0000    0.9994    0.9976    0.9945    0.9903    0.9848    0.9781    0.9703    0.9613    0.9511    0.9397...
    0.9272    0.9135    0.8988    0.8829    0.8660    0.8480    0.8290    0.8090    0.7880    0.7660    0.7431    0.7193    0.6947    0.6691...
    0.6428    0.6157    0.5878    0.5592    0.5299    0.5000    0.4695    0.4384    0.4067    0.3746    0.3420    0.3090    0.2756    0.2419...
    0.2079    0.1736    0.1392    0.1045    0.0698    0.0349    0.0000   -0.0349   -0.0698   -0.1045   -0.1392   -0.1736   -0.2079   -0.2419...
    -0.2756   -0.3090   -0.3420   -0.3746   -0.4067   -0.4384   -0.4695   -0.5000   -0.5299   -0.5592   -0.5878   -0.6157   -0.6428   -0.6691...
    -0.6947   -0.7193   -0.7431   -0.7660   -0.7880   -0.8090   -0.8290   -0.8480   -0.8660   -0.8829   -0.8988   -0.9135   -0.9272   -0.9397...
    -0.9511   -0.9613   -0.9703   -0.9781   -0.9848   -0.9903   -0.9945   -0.9976   -0.9994   -1.0000   -0.9994   -0.9976   -0.9945   -0.9903...
    -0.9848   -0.9781   -0.9703   -0.9613   -0.9511   -0.9397   -0.9272   -0.9135   -0.8988   -0.8829   -0.8660   -0.8480   -0.8290   -0.8090...
    -0.7880   -0.7660   -0.7431   -0.7193   -0.6947   -0.6691   -0.6428   -0.6157   -0.5878   -0.5592   -0.5299   -0.5000   -0.4695   -0.4384...
    -0.4067   -0.3746   -0.3420   -0.3090   -0.2756   -0.2419   -0.2079   -0.1736   -0.1392   -0.1045   -0.0698   -0.0349   -0.0000];

myInd = theta_to_ind_for_Cos_and_Sin(wrapped_input_angle, length(SinTable));
mySin = SinTable(myInd);
end


function myInd = theta_to_ind_for_Cos_and_Sin(wrapped_input_angle, myLength)
myInd = round(1 + (wrapped_input_angle/ (2*pi))* myLength);
if myInd > myLength
    myInd = myLength;
end
end

function wrapped_input_angle = mywrapTo2Pi(myInput, app_level);
% This is copied from Matlab's
positiveInput = (myInput > 0);

if ~isempty(strfind(app_level, 'B'))
    myPI = fi(pi, 1, 18, 17);
elseif ~isempty(strfind(app_level, 'A'))
    myPI = fi(pi, 1, 27, 12);
else
    myPI = pi;
end
myInput = mod(myInput, 2*myPI);
myInput((myInput == 0) & positiveInput) = 2*pi;
wrapped_input_angle = myInput;

end