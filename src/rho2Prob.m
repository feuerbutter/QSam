function prob = rho2Prob(rho,POM)
    % 
    % This calculates the observation probabilities corresponding to the 
    % measurement operators, POM, for the states, rho.
    % 
    % Input
    % --------------------------------------------------------------------------
    % rho : array of complex double
    %   density matrices
    % POM : array of complex double
    %   measurement operators
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % prob : array of real double
    %   observation probabilities
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
    
    % Flatten rho into columns
    n_rho = size(rho,3);
    rho_c = reshape(rho,m^2,n_rho);

    % calculate prob
    prob = real(POM_r*rho_c);
    