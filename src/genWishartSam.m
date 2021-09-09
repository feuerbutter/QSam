function [rhos,Sigma] = genWishartSam(N,rho_peak,n_w)
    % 
    % This function generates a Wishart sample that is peaked at rho_peak 
    % and has N sample points using Psi matrices with n_w columns.
    % Input
    % --------------------------------------------------------------------------
    % N : int
    %   # of sample points generated
    % rho_peak : m*m array of complex
    %   peak of the Wishart distribuion
    % n_w : int
    %   # of columns of the Psi matrices in constructing the density
    %   matrices
    % 
    % Output
    % --------------------------------------------------------------------------
    % rhos : 3d array of comnplex double
    %   sample points in the state space
    % 
    % 
    % 
    
    m = size(rho_peak,1);
    Sigma = calSigma(n_w,rho_peak);
    A = mpower(Sigma,1/2);

    Psi_real = randn([m,n_w,N]); 
    Psi_imag = randn([m,n_w,N]); 
    Psi = Psi_real + 1j * Psi_imag;

    rhos = zeros(m,m,N);
    
    for n_dx = 1 : N
       rhotemp = Psi(:,:,n_dx) * (Psi(:,:,n_dx))' ; 
       rhotemp =  A * rhotemp * A;
       rhotemp = rhotemp / trace(rhotemp); 
       rhos(:,:,n_dx) = rhotemp;
    end
    