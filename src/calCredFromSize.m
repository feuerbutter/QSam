function cred_from_size = calCredFromSize(size_lambda,log_lambda_step,min_log_lambda,max_log_lambda)
    % 
    % This function calculates the credibility theoretically from the size 
    % according to Eqn.(56), 
    % given that the size can be approximated to a high precision
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
    
    % calculate the log of lambdas
    log_lambdas = min_log_lambda:log_lambda_step:max_log_lambda;
    n_lambda = length(log_lambdas);

    % the lambda values
    lambdas = exp(log_lambdas);
    % normalisation factor of Eqn.(56)
    nors = sum(size_lambda.*lambdas');
    
    % initialize the credibility variable
    cred_from_size = zeros(n_lambda,1);
    cred_from_size(1) = 1;

    % calculating cred according to Eqn.(56), 
    for l_dx = 2:n_lambda-1
        cred_from_size(l_dx) = (size_lambda(l_dx)*lambdas(l_dx)+sum(size_lambda(l_dx:end).*lambdas(l_dx:end)')*log_lambda_step)/(nors*log_lambda_step);
    end

    cred_from_size(n_lambda) = 0;
    