% The main run script for QSam
tic
%------- Properties of the reference sample
% Total number of reference samples to produce
N_total = 1e7;

kappa_w = 0.; % percentage of Wishart points
n_w = 5; % number of columns for the Wishart sample

kappa_s = 0.9; %percentage of linearly-shifted points
n_s = 10;

lambso = 0.8; %the 'distance' between the origin and the ref peak before linear shift, 
% 0 meaning the ref peak is at the origin, 1 at MLE
lambs = 0.95; %to which place the ref sample is shifted, 1 being MLE, 0 being origin

lambp = 0.95; %the 'distance' between the origin and the Wishart ref peak, 
% 0 meaning the Wishart ref peak is at the origin, 1 at MLE

% Target posterior counts
pop = 10 * [10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6];

% corpacc contains the output target sample and ar is the acceptance rate.
[corpacc,ar] = genTarSam(N_total,pop,kappa_w,n_w,lambp,kappa_s,n_s,lambso,lambs);
toc
