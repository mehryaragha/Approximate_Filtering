% Written by: Mehryar Emambakhsh
% Email: mehryar_emam@yahoo.com
% Date: 27 Sep 2017
% Paper:
% P. Garcia, M. Emambakhsh, A. Wallace, “Learning to Approximate Computing at Run-time,”
% in IET 3rd International Conference on Intelligent Signal Processing (ISP
% 2017), 2017, to appear.

% This is a Demo file, simulating Figure 4 of the paper above.

% Exact computation
Demo_Exact_Computation;

% Kullback–Leibler (KL) threshold 0.5
my_KL_thresh = 0.5; Dynamic_approximation(my_KL_thresh);

% KL threshold 1
my_KL_thresh = 1;
Dynamic_approximation(my_KL_thresh);

% KL threshold 2
my_KL_thresh = 2;
Dynamic_approximation(my_KL_thresh);

% KL threshold 4
my_KL_thresh = 4;
Dynamic_approximation(my_KL_thresh);
