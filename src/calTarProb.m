function [log_tar_prob,prob_points] = calTarProb(rhos,pop,pom)
    % 
    % This calculates the log of target probability for the case of posterior
    % 
    % 
    % Input
    % --------------------------------------------------------------------------
    % rhos : 3d array of comnplex double
    %   sample points in the state space
    % pop : array of integers
    %   the power of probabilities for the posterior
    % pom : array of complex double
    %   probability operator measurement
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % log_tar_prob : 1d array of real double
    %   log of target probabilities of the sample
    % prob_points : 2d array of complex double 
    %   sample points in the probability space
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see genTarSam.m
    % 
    % 
    % 
    
    % extract the number of states and the dimension of the system
    N = size(rhos,3);
    m = size(rhos,1);
    
    % calculate the probabilities, where each sample probability is
    % returned in a column, hence, prob_points is a m^2 by N matrix
    prob_points = rho2Prob(rhos,pom);
    
    % pop is input as a row and popl is N copies of transpose of pop
    popl = repmat(pop',1,N);
    
    
    % This element-wise multiplication calculates the log of (the power of 
    % the probabilities)
    log_tar_prob = log(prob_points) .* popl;
    
    % The sum gives the log of posterior
    log_tar_prob = sum(log_tar_prob);
