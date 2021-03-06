function [rhos,num_phys,Sigma_w,Sigma_s] = genRefSam(n_w,n_s,N_w,N_s,N_total,m,rho_peak,rho_peak_shift_o,rho_peak_shift)
    % 
    % This would generate a reference sample with Nt total points, in which 
    % there are Nw Wishart sample points peaked at rho_peak, Ns linearly shifted 
    % Wishart sample points, which is generated by linearly shifting a Wishart 
    % distribution peaked at rho_peak_shift_o to rho_peak_shift and (Nt-Nw-Ns) uniformly generated
    % background points.
    % 
    % Input
    % --------------------------------------------------------------------------
    % n_w : int
    %   # of columns of the Psi matrix used in construction of the Wishart sample
    % n_s : int
    %   # of columns of the Psi matrix used in construction of the shifted sample
    % N_w : int
    %   # of Wishart sample points
    % N_s : int
    %   # of shifted sample points
    % N_total : int
    %   # of sample points generated
    % m : int
    %   dimension of the system, e.g. 2^k for k qubits
    % rho_peak : array of complex double
    %   peak of the Wishart sample
    % rho_peak_shift_o : array of complex double
    %   peak of the shifted sample before shifting
    % rho_peak_shift : array of complex double
    %   the final peak of the shifted sample
    % 
    % Output
    % --------------------------------------------------------------------------
    % rhos : 3d array of comnplex double
    %   sample points in the state space
    % num_phys : int
    %   # of physical states in the reference sample, only different from the total number of points 
    %   in the case of a shifted sample
    % Sigma_w : array of double
    %   the covariance matrix required to produce the peaked Wishart sample
    % Sigma_s : array of double
    %   the covariance matrix required to produce the shifted Wishart sample
    % 
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see genTarSam.m
    % 
    % 
    % 
    
    N_uni = N_total-N_w-N_s;

    rhos = zeros(m,m,N_total);
    
    %-----wishart
    [rho_w,Sigma_w] = genWishartSam(N_w,rho_peak,n_w);
    
    %-----adding a uniform background
    rho_u = genUniSam(N_uni,m);
    
    %-----shifted Wishart
    [rho_s,num_phys,Sigma_s] = genShiftedSam(N_s,rho_peak_shift_o,rho_peak_shift,n_s);
    
    rhos(:,:,1:N_w) = rho_w;
    rhos(:,:,N_w+1:N_w+N_uni) = rho_u;
    rhos(:,:,N_w+N_uni+1:N_total) = rho_s;
    
    