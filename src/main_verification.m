% Script for verification by plotting size and credibility
tic
rng('shuffle');

%% Parameters Setting

%--- General
Nt = 1e5; %number of samples in a run, better not go beyond 5e5, otherwise the computer might hang
Nrun = 2; %number of runs in a trial
GetCorp = 0; %choice to save the sample, 1 for both the sample and the 
             %acceptance rate, 0 for only the acceptance rate

NT = Nrun*Nt; %total number of reference sample points

nq = 2; % # of qubits
m = 2^nq; % for qubits only, if qutrit, modified d and pom acoordingly
np = m^2;
pom = buildNTetraPOM(nq); % the default measurement is multiple tetrahedron
% pop = [2,2,2,2]; % target posterior
% pop = [10,20,25,45];
pop = 10 * [10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6];

%  pow = 2;
%  pop = ones(1,np)*pow; % column of np probabilities

% pop = [36,13,64,71,14,16,7,15,60,10,84,63,64,9,55,71,8,12,10,16,16,48,67,62,9,64,75,63,10,74,60,73,65,14,62,66,9,57,76,53,82,78,128,22,61,44,25,27,56,12,52,66,14,76,56,78,45,47,22,27,66,68,25,102];


%--- Reference Distribution
kappa_w = 0.; %percentage of Wishart-shifted points
n_w = 80; 

kappa_s = 0.6; %percentage of linearly-shifted points
n_s = 5;

lambso = 0.75; %the 'distance' between the origin and the ref peak before linear shift, 
% 0 meaning the ref peak is at the origin, 1 at MLE
lambs = 0.9; %to which place the ref sample is shifted, 1 being MLE, 0 being origin

lambp = 0.95; %the 'distance' between the origin and the Wishart ref peak, 
% 0 meaning the Wishart ref peak is at the origin, 1 at MLE


%% Sampling Procedure

%--- MLE
np = m^2;
ref_freq = pop/sum(pop);

rho_peak = prob2Rho(ref_freq,pom);
v = eig(rho_peak);

if sum(v<0)
    fprintf('linear inversion not physical, mle running \n');
    rhomle = qse_apg(pom,ref_freq');
    fprintf('mle done \n');
    rho_peak = lambp*rhomle+(1-lambp)*eye(m)/m;   
    ppeak = rho2Prob(rhomle,pom)';
    %-- o. peak for the shifted ref
    rho_peak_shift_o = rhomle*lambso+(1-lambso)*eye(m)/m;
    rho_peak_shift = rhomle*lambs+(1-lambs)*eye(m)/m;
else
    rho_peak_shift_o = rho_peak*lambso+(1-lambso)*eye(m)/m;
    rho_peak_shift = rho_peak*lambs+(1-lambs)*eye(m)/m; 
    ppeak = rho2Prob(rho_peak,pom)';
end  
    
logLpeak = sum(pop.*log(ppeak));
%-----Plotting
log_lambda_step = 0.01;
min_log_lambda = -20;
max_log_lambda = 0;

% generation of a large uniform sample to calculate size
N_size_uni = 1e7;
rho_u = genUniSam(N_size_uni,m);
prob_points_uniform = rho2Prob(rho_u,pom);
size_lambda = calCred(pop,prob_points_uniform,logLpeak,log_lambda_step,min_log_lambda,max_log_lambda);
cred_from_size = calCredFromSize(size_lambda,log_lambda_step,min_log_lambda,max_log_lambda);
cred_lambda = calCred(pop,prob_points_accepted,logLpeak,log_lambda_step,min_log_lambda,max_log_lambda);
cred_py = calCred(pop,prob_points_phys',logLpeak,log_lambda_step,min_log_lambda,max_log_lambda);

%------Plotting for verification
plotCS(cred_from_size,log_lambda_step,min_log_lambda,max_log_lambda)
hold on
plotCS(cred_lambda,log_lambda_step,min_log_lambda,max_log_lambda)
plotCS(cred_py,log_lambda_step,min_log_lambda,max_log_lambda)
