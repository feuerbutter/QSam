function cred_lambda = calCred(pop,prob_points,logLpeak,log_lambda_step,min_log_lambda,max_log_lambda)
    % 
    % This calculates the credibility curve of the sample points in the probability space, 
    % given in 'prob_points' for range of the log of lambda from min_log_lambda to 
    % max_log_lambda in steps of log_lambda_step
    % 
    % Input
    % --------------------------------------------------------------------------
    % pop : array of integers
    %   the power of probabilities for the posterior
    % prob_points : 3d array of complex double
    %   sample points in the probability space, e.g. a 2-qubit sample of 1000 points would
    %   have the dimension of 4x4x1000
    % logLpeak : double
    %   the log of the maximum likelihood in the physical space 
    % log_lambda_step : double
    %   the x-axis interval between the points forming the curve
    % min_log_lambda : double
    %   the x-axis starting point of the curve
    % max_log_lambda : double 
    %   the x-axis ending point of the curve, usually set to 0 = log(1)
    % 
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % cred_lambda : 1d array of double
    %   the credibility values for all the lambdas
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see main_verification.m
    % 
    % 

    % pop is given in a row
    % Each point of prob_points is given in a column
    % extract the number of points in the sample
    N = size(prob_points,2);

    % calculate the log of lambdas
    log_lambdas = min_log_lambda:log_lambda_step:max_log_lambda;
    n_lambda = length(log_lambdas);

    % credibility at each lambda value
    cred_lambda = zeros(n_lambda,1);
    
    % pop is input as a row and popl is N copies of transpose of pop
    popl = repmat(pop',1,N);

    % This element-wise multiplication calculates the log of (the power of 
    % the probabilities), and thus the sum gives the log of the likelihood as 
    % given in Eqn.(43) 
    log_likelihood = sum(log(prob_points).*popl);
    
    % This calculates the credibility according to Eqn.(54) by counting how
    % many points with a likelihood larger than lambda * (maximum likelihood)
    for l_dx = 1:n_lambda
        cred_lambda(l_dx) = sum(log_likelihood>(log_lambdas(l_dx)+logLpeak))/N;
    end
    