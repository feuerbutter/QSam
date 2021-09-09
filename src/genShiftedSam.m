function [rho,num_phys,Sigma] = genShiftedSam(N,rho_peak_shift_o,rho_peak_shift,n_s)
    % 
    % This function generates a Wishart sample that is peaked at rho_peak 
    % and has N sample points using Psi matrices with n_w columns.
    % Input
    % --------------------------------------------------------------------------
    % N : int
    %   # of sample points generated
    % rho_peak_shift_o : array of complex double
    %   peak of the shifted sample before shifting
    % rho_peak_shift : array of complex double
    %   the final peak of the shifted sample
    % n_s : int
    %   # of columns of the Psi matrices in constructing the density
    %   matrices
    % 
    % Output
    % --------------------------------------------------------------------------
    % rhos : 3d array of comnplex double
    %   sample points in the state space
    % num_phys : int
    %   # of physical states in the shifted reference sample
    % Sigma : array of double
    %   the covariance matrix required to produce the shifted Wishart sample
    % 
    % 
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see genRefSam.m
    % 
    % 
    % 
    
    m = size(rho_peak_shift_o,1);
    delrho = rho_peak_shift-rho_peak_shift_o;

    Sigma = calSigma(n_s,rho_peak_shift_o);
    A = mpower(Sigma,1/2);

    Psi_real = randn([m,n_s,N]); 
    Psi_imag = randn([m,n_s,N]); 
    Psi = Psi_real + 1j * Psi_imag;

    rho = zeros(m,m,N);

    num_phys = 0; % # of physical states after shifting

    for n_dx = 1 : N
        rhotemp = Psi(:,:,n_dx) * (Psi(:,:,n_dx))' ; 
        rhotemp =  A * rhotemp * A;
        rhotemp = rhotemp / trace(rhotemp); 
        rhotemp = rhotemp + delrho;
        if min(real(eig(rhotemp)))<0
            continue;
        end
        num_phys = num_phys + 1;
        rho(:,:,num_phys) = rhotemp;
    end
    
