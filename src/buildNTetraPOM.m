function n_tetra_pom = buildNTetraPOM(n_qb)
    % 
    % Build the tetrahedron POM.
    % 
    % Input
    % --------------------------------------------------------------------------
    % n_qb : int
    %   # of qubits
    % 
    % Output
    % --------------------------------------------------------------------------
    % pom : array of real double
    %   tetrahedron pom for n_qb qubits
    % 
    % How to call
    % --------------------------------------------------------------------------
    % n_qb = 2
    % n_tetra_pom = buildNTetraPOM(n_qb) % the tensor product of 2 tetrahedron poms to form a pom in the 2 qubits space
    % 


    % pauli matrices
    sigma = zeros(2,2,3);
    sigma(:,:,1) = [0,1;1,0];
    sigma(:,:,2) = [0,-1i;1i,0];
    sigma(:,:,3) = [1,0;0,-1];

    % build the fully symmetric tetrahedron POM for 1 qubit
    pom(:,:,1) = eye(2)/4 + 1/4/sqrt(3) * (sigma(:,:,1)-sigma(:,:,2)-sigma(:,:,3));
    pom(:,:,2) = eye(2)/4 + 1/4/sqrt(3) * (-sigma(:,:,1)+sigma(:,:,2)-sigma(:,:,3));
    pom(:,:,3) = eye(2)/4+1/4/sqrt(3) * (-sigma(:,:,1)-sigma(:,:,2)+sigma(:,:,3));
    pom(:,:,4) = eye(2)/4+1/4/sqrt(3) * (sigma(:,:,1)+sigma(:,:,2)+sigma(:,:,3));

    po = pom;
            
    if n_qb > 1    
        for q_dx = 1:n_qb-1
            tempn = 2^(2*q_dx);
            tempom = zeros(2^(q_dx+1),2^(q_dx+1),tempn*4);
            for j_dx = 1:4
                for t_dx = 1:tempn
                    tempom(:,:,tempn*(j_dx-1)+t_dx) = kron(pom(:,:,t_dx),po(:,:,j_dx));
                end
            end
            pom = tempom;
        end
    else
    
    end
    n_tetra_pom = pom;

