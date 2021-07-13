function Sigma = calSigma(n,rho_peak)
    % 
    % This calculates the covariance matrix, Sigma, according to Eqn.(25)
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

    m = size(rho_peak,1);

    % according to Eqn.(25)
    Sigma = eye(m) / (inv(rho_peak)+eye(m)*m^2/(n-m));
    