function rho = genUniSam(N,m)
    % 
    % This generates a uniform sample in m dimension with N points.
    % 
    % Input
    % --------------------------------------------------------------------------
    % N : int
    %   # of sample points generated
    % m : int
    %   dimension of the system, e.g. 2^k for k qubits
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % 
    % 
    % 
    % 

    Psi_real = randn([m,m,N]); 
    Psi_imag = randn([m,m,N]); 
    Psi = Psi_real + 1j * Psi_imag;
    
    rho = zeros(m,m,N);

    for n_dx = 1 : N
       rhotemp = Psi(:,:,n_dx)*(Psi(:,:,n_dx))' ; 
       rho(:,:,n_dx) = rhotemp/trace(rhotemp); 
    end  
    
