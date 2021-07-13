function [rho,num_phys,Sigma] = genShiftedSam(N,rho_peak_shift_o,rho_peak_shift,n_s)
    % 
    % This generates a Wishart sample originally peaked at rho_peak_shift_o and
    % is shifted to rho_peak_shift.
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
    
