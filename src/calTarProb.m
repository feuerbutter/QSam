function [log_tar_prob,corp] = calTarProb(rhos,pop,pom)
    % 
    % This calculates the log of target probability for the case of posterior
    % 
    % 
    % Input
    % --------------------------------------------------------------------------
    % rhos : 3d array of complex
    %   the samples
    % pop : row array of real
    %   the posterior
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % 
    % 
    % 
    % 
    
    % extract the number of states and the dimension of the system
    N = size(rhos,3);
    m = size(rhos,1);
    
    % calculate the probabilities, where each sample probability is
    % returned in a column, hence, corp is a m^2 by N matrix
    corp = rho2Prob(rhos,pom);
    
    % pop is input as a row and popl is N copies of transpose of pop
    popl = repmat(pop',1,N);
    
    
    % This element-wise multiplication calculates the log of (the power of 
    % the probabilities)
    log_tar_prob = log(corp) .* popl;
    
    % The sum gives the log of posterior
    log_tar_prob = sum(log_tar_prob);
