function log_w_prob = calWishartProb(rho,Ncol,m,inv_Sigma)
    % 
    % This calculates the logarithm of the part of the Wishart probability 
    % that has a dependence on rho according to Eqn.(16)
    % 
    % Input
    % --------------------------------------------------------------------------
    % rho : 2d array of complex double
    %   1 single state matrix, as the det function does not take in a stack of matrices
    % n_w : int
    %   # of columns of the Psi matrix used in construction of the Wishart sample
    % m : int
    %   dimension of the system, e.g. 2^k for k qubits
    % inv_Sigma : array of double
    %   inverse of the covariance matrix required to produce the peaked Wishart sample
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % log_w_prob : real
    %   log of the part of the Wishart probability that depends on rho
    % 
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see calRefProb.m
    % 
    
    log_w_prob = real((Ncol-m)*log(det(rho))-(Ncol*m)*log(trace(rho*inv_Sigma)));
    