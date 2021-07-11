function pom = buildTetraPOM()
    % 
    % Build the tetrahedron POM.
    % 
    % Input
    % --------------------------------------------------------------------------
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % pom : array of real double
    %   tetrahedron pom
    % 
    % 
    % 

    % pauli matrices
    sigma = zeros(2,2,3);
    sigma(:,:,1) = [0,1;1,0];
    sigma(:,:,2) = [0,-1i;1i,0];
    sigma(:,:,3) = [1,0;0,-1];

    % build the fully symmetric tetrahedron POM
    pom(:,:,1) = eye(2)/4 + 1/4/sqrt(3) * (sigma(:,:,1)-sigma(:,:,2)-sigma(:,:,3));
    pom(:,:,2) = eye(2)/4 + 1/4/sqrt(3) * (-sigma(:,:,1)+sigma(:,:,2)-sigma(:,:,3));
    pom(:,:,3) = eye(2)/4+1/4/sqrt(3) * (-sigma(:,:,1)-sigma(:,:,2)+sigma(:,:,3));
    pom(:,:,4) = eye(2)/4+1/4/sqrt(3) * (sigma(:,:,1)+sigma(:,:,2)+sigma(:,:,3));
    