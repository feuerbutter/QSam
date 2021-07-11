function cred_lambda = calCred(pop,corp,logLpeak,log_lambda_step,min_log_lambda,max_log_lambda)
    % 
    % This calculates the credibility in the samples in the probability space, 
    % given in 'corp' for range of the log of lambda from min_log_lambda to 
    % max_log_lambda in steps of log_lambda_step
    % 
    % Input
    % --------------------------------------------------------------------------
    % n : int
    %   # of columns of the A matrix used in construction of the Wishart sample
    % rho_peak : array of complex double
    %   peak of the Wishart sample
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % 
    % 
    % 
    % 

    % pop is given in a row
    % Each point of Corp is given in a column
    % extract the number of 
    N = size(corp,2);

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
    log_likelihood = sum(log(corp).*popl);
    
    % This calculates the credibility according to Eqn.(54) by counting how
    % many points with a likelihood larger than lambda * (maximum likelihood)
    for l_dx = 1:n_lambda
        cred_lambda(l_dx) = sum(log_likelihood>(log_lambdas(l_dx)+logLpeak))/N;
    end
    