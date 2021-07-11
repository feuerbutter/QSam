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
pop = [10; 4; 6; 4; 7; 6; 5; 6; 5; 6; 10; 6; 5; 6; 8; 6]';

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

rhop = prob2Rho(ref_freq,pom);
v = eig(rhop);

if sum(v<0)
    fprintf('linear inversion not physical, mle running \n');
    rhomle = qse_apg(pom,ref_freq');
    fprintf('mle done \n');
    rhop = lambp*rhomle+(1-lambp)*eye(m)/m;   
    ppeak = rho2Prob(rhomle,pom)';
    %-- o. peak for the shifted ref
    rhopso = rhomle*lambso+(1-lambso)*eye(m)/m;
    rhops = rhomle*lambs+(1-lambs)*eye(m)/m;
else
    rhopso = rhop*lambso+(1-lambso)*eye(m)/m;
    rhops = rhop*lambs+(1-lambs)*eye(m)/m; 
    ppeak = rho2Prob(rhop,pom)';
end  
    
logLpeak = sum(pop.*log(ppeak));
    

%% --- Resampling

Nw = ceil(Nt*kappa_w);
Ns = floor(Nt*kappa_s);
Nb = Nt-Nw-Ns;
kappa = 1 - kappa_w - kappa_s;

lalall = zeros(Nt,Nrun);
corpall = zeros(m^2,Nt,Nrun);
pNt = zeros(Nrun,1);
    
for nrunk = 1:Nrun
%     for nrunk = 1:Nrun

    %-- obtaining the ref sam
    [rhos,num_phys,Sigma_w,Sigma_s] = genRefSam(n_w,n_s,Nw,Ns,Nt,m,rhop,rhopso,rhops);
    %corp/rhos returns the sample in the p/rho-space
    %CVM are the convariance matrices
    %indunph records the index of the points that could only be from 
    %the Wishart but not the linear shift part so that the resampling
    %ratio can be calculated correctly

    %-- calculation of the reference probability of the points
    ref_prob = calRefProb(rhos,kappa,kappa_w,n_w,Sigma_w,n_s,Sigma_s,rhopso,rhops);
%     ref_prob = calRefProb(rhos,kappa,kappa_w,n_w,Sigma_w);

    %-- calculation of the target probability of the points
    [lfx,corp] = calTarProb(rhos,pop,pom);

    %-- the acceptance ratio in log is "lal"
    lrx = log(ref_prob);
    lal = lfx - lrx;

    lalall(:,nrunk) = lal;
    corpall(:,:,nrunk) = corp;
end

toc



lallin = lalall(:);
[lalcombdam,inal] = max(lallin);
corpalllin = reshape(corpall,[np,NT]);
ldraw = log(rand(size(lallin)));
acc = ldraw < (lallin-lalcombdam);
corpacc = corpalllin(:,acc);
ar = size(corpacc,2)/NT;
toc

%-----Plotting
log_lambda_step = 0.01;
min_log_lambda = -20;
max_log_lambda = 0;

N_size_uni = 1e7;
corpu = genUniSam(N_size_uni,m);
corpu = rho2Prob(corpu,pom);
size_lambda = calCred(pop,corpu,logLpeak,log_lambda_step,min_log_lambda,max_log_lambda);
cred_from_size = calCredFromSize(size_lambda,log_lambda_step,min_log_lambda,max_log_lambda);
cred_lambda = calCred(pop,corpacc,logLpeak,log_lambda_step,min_log_lambda,max_log_lambda);

%------Plotting for verification
plotCS(cred_from_size,log_lambda_step,min_log_lambda,max_log_lambda)
hold on
plotCS(cred_lambda,log_lambda_step,min_log_lambda,max_log_lambda)

