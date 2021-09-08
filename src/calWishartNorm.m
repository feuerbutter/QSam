function log_w_norm = calWishartNorm(n_w,m,inv_Sigma)
    % 
    % This calculates the logarithm of the part of the Wishart probability 
    % that is independent of rho according to Eqn.(16), which is essentially
    % the normalization of the distribution
    % 
    % Input
    % --------------------------------------------------------------------------
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
    % log_w_norm : real double
    %   log of the normalization
    %   
    % How to call
    % --------------------------------------------------------------------------
    % see calRefProb.m
    % 
    
    gl = gammaln(n_w+1-(1:m));
    
    log_w_norm = gammaln(n_w*m) - m * (m-1) / 2 * log(pi) - sum(gl) + n_w * log(real(det(inv_Sigma)));
