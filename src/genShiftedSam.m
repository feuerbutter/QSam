function [rho,num_phys,Sigma] = genShiftedSam(N,rho_peak_shift_o,rho_peak_shift,n_s)
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
    % 
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
    
