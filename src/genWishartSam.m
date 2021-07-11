function [rho,Sigma] = genWishartSam(N,rhop,n_w)
    % 
    % This generates a uniform sample in m dimension with N points. If one 
    % has a specific POM in mind, one can set prob_opt to be true and have the
    % points in the probability space returned as well.
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
    
    m = size(rhop,1);
    Sigma = calSigma(n_w,rhop);
    kra = mpower(Sigma,1/2);

    Psi_real = randn([m,n_w,N]); 
    Psi_imag = randn([m,n_w,N]); 
    Psi = Psi_real + 1j * Psi_imag;

    rho = zeros(m,m,N);
    
    for n_dx = 1 : N
       rhotemp = Psi(:,:,n_dx) * (Psi(:,:,n_dx))' ; 
       rhotemp =  kra * rhotemp * kra;
       rhotemp = rhotemp / trace(rhotemp); 
       rho(:,:,n_dx) = rhotemp;
    end
    