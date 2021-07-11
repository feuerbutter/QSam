function log_w_norm = calWishartNorm(Ncol,m,CVMin)
    % 
    % This calculates the logarithm of the part of the Wishart probability 
    % that is independent of rho according to Eqn.(16), which is essentially
    % the normalization of the distribution
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
    
    gl = gammaln(Ncol+1-(1:m));
    
    log_w_norm = gammaln(Ncol*m) - m * (m-1) / 2 * log(pi) - sum(gl) + Ncol * log(real(det(CVMin)));
