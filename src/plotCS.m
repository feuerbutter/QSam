function plotCS(size_or_cred,log_lambda_step,min_log_lambda,max_log_lambda)
    % 
    % This is the utility function to draw a semilogx plot of size/credibility
    % against lambda
    % 
    % Input
    % --------------------------------------------------------------------------
    % n : int
    %   # of columns of the Psi matrix used in construction of the Wishart sample
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
    
    log_lambdas = min_log_lambda:log_lambda_step:max_log_lambda;
    plot(log_lambdas,size_or_cred);