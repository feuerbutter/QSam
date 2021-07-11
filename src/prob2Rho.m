function rho = prob2Rho(prob,POM)
    % 
    % This calculates the density matrix (might not be physical) corresponding
    % to the observation probability or relative frequency, prob, with the
    % measurement, POM.
    % 
    % Input
    % --------------------------------------------------------------------------
    % prob : array of real double
    %   observation probability or relative frequency
    % POM : array of complex double
    %   measurement operators
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % rho : array of complex double
    %   density matrices
    % 
    % 
    % 

    % Flatten POM
    n_POM = size(POM,3);
    m = size(POM,2);
    
    % transpose each POM
    POM_trans = permute(POM,[2 1 3]);
    
    % flatten each POM into a row
    POM_c = reshape(POM_trans,m^2,n_POM);
    POM_r = POM_c.';

    rho_c = POM_r\prob.';
    rho = reshape(rho_c,m,m);
