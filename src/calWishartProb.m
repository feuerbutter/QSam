function log_w_prob = calWishartProb(rho,Ncol,m,inv_Sigma)
    % 
    % This calculates the logarithm of the part of the Wishart probability 
    % that has a dependence on rho according to Eqn.(16)
    % 
    % Input
    % --------------------------------------------------------------------------
    % N : int
    %   # of sample points generated
    % m : int
    %   dimension of the system, e.g. 2^k for k qubits
    % varargin : arrays of complex
    %   POM is supplied as an optional variable
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % 
    % 
    % 
    % 
    
    log_w_prob = real((Ncol-m)*log(det(rho))-(Ncol*m)*log(trace(rho*inv_Sigma)));
    