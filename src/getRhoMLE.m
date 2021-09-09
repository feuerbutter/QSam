function rho_mle = getRhoMLE(pop,pom)
    % 
    % Find the physical maximum likehood state for a given posterior, given by 
    % pop, and an informationally complete set of measurement, pom.
    % 
    % 
    % Input
    % --------------------------------------------------------------------------
    % pop : array of integers
    %   the power of probabilities for the posterior
    % pom : array of complex double
    %   probability operator measurement
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % rho_mle : 2d array of complex double
    %   the physical maximum likehood state 
    % 
    % 
    % 

    ref_freq = pop/sum(pop);

    % Linearly invert the probabilities to get the state
    rho_lin_inv = prob2Rho(ref_freq,pom);

    % Calculate the eigenvalues to check physicality
    v = eig(rho_lin_inv);

    if sum(v<0)
        % rho_lin_inv not physical
        fprintf('linear inversion not physical, mle running \n');
        rho_mle = qse_apg(pom,ref_freq');
        fprintf('mle done \n');
    else
        rho_mle = rho_lin_inv;
    end
