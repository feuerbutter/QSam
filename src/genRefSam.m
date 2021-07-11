function [rhos,num_phys,Sigma_w,Sigma_s] = genRefSam(n_w,n_s,N_w,N_s,N_total,m,rhop,rhopso,rhops)
    % 
    % This would generate a reference sample with Nt total points, in which 
    % there are Nw Wishart sample points peaked at rhop, Ns linearly shifted 
    % Wishart sample points, which is generated by linearly shifting a Wishart 
    % distribution peaked at rhopso to rhops and (Nt-Nw-Ns) uniformly generated
    % background points.
    % 
    % Input
    % --------------------------------------------------------------------------
    % n_w : int
    %   # of columns of the A matrix used in construction of the Wishart sample
    % n_s : int
    %   # of columns of the A matrix used in construction of the shifted sample
    % Nw : int
    %   # of Wishart sample points
    % Ns : int
    %   # of shifted sample points
    % Nt : int
    %   # of total sample points
    % d : int
    %   dimension of the system, e.g. 2^k for k qubits
    % rhop : array of complex double
    %   peak of the Wishart sample
    % rhops : array of complex double
    %   the final peak of the shifted sample
    % rhopso : array of complex double
    %   peak of the shifted sample before shifting
    % TransM : ?????
    %   POVM ????
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % 
    % 
    % 
    % 

    
    N_uni = N_total-N_w-N_s;

    rhos = zeros(m,m,N_total);
    
    %-----wishart
    [rho_w,Sigma_w] = genWishartSam(N_w,rhop,n_w);
    
    %-----adding a uniform background
    rho_u = genUniSam(N_uni,m);
    
    %-----shifted Wishart
    [rho_s,num_phys,Sigma_s] = genShiftedSam(N_s,rhopso,rhops,n_s);
    
    rhos(:,:,1:N_w) = rho_w;
    rhos(:,:,N_w+1:N_w+N_uni) = rho_u;
    rhos(:,:,N_w+N_uni+1:N_total) = rho_s;
    
    