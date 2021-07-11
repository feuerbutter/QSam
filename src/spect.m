%%%%%%
function [A,JacDet,u]=spect(q,d,pom)

    % d=4, for example
    nt=d*(d+1)/2-1; % number of theta parameters in HMC
    nf=d*(d-1)/2; % number of phi parameters in HMC
    num=nt+nf; % total number of indpendent parameters
    
    t=q(1:nt)';
    f=q(nt+1:num)';
    
    sint=sin(t); cost=cos(t);
    tant=tan(t); % cott=cot(t);
    expf=exp(-1i*f);
    
    % product of sin(t)
    % sp(1)=1; sp(k)=sint(1)*sint(2)*...*sint(k-1) for k>1
    sp=[1 cumprod(sint)];
    % Create x, which are |Ajk| and lie on a sphere
    x=sp.*[cost 1];
    
    % Create A matrix, A is upper triangular
    % rho=A^\dagger * A
    % Entries are Ajk=0 if j<k;
    % A11=x(1), A12=exp(1i f(1)) x(2), A13=exp(1i f(2)) x(3), etc.; 
    % Diagonals, Aii, are real.
    A=zeros(d); % construct matrix A with x (theta) and phi
    xt=zeros(nt+1,2); % j,k indices for theta
    ind1=1;
    for j=1:d
        for k=j:d
            xt(ind1,:)=[j,k];
            A(j,k)=x(ind1);
            ind1=ind1+1;
        end
    end
    
    ft=zeros(nf,2); % j,k indices for phi
    ind2=1;
    for j=1:(d-1)
        for k=(j+1):d
            ft(ind2,:)=[j,k];
            A(j,k)=A(j,k)*expf(ind2);
            ind2=ind2+1;
        end
    end
    
    cA=A'; % complex transpose of A
    
    %
    dAdt=zeros(d,d,nt); % derivative of A wrt theta: dA/dt(m)
    dAdtdf=zeros(d,d,nt,nf); % derivative of A wrt theta then phi
    for m=1:nt % theta
        ind=xt(m,:);
        dAdt(ind(1),ind(2),m)=-A(ind(1),ind(2))*tant(m);
        
        for k=(m+1):(nt+1)
            ind=xt(k,:);
            dAdt(ind(1),ind(2),m)=A(ind(1),ind(2))/tant(m);
        end
        
        for j=1:nf
            ind=ft(j,:);
            dAdtdf(ind(1),ind(2),m,j)=-1i*dAdt(ind(1),ind(2),m);
        end
    end
    
    % 2nd order derivative wrt theta
    ddAddt=zeros(d,d,nt,nt);
    for i=1:nt
        for k=i:nt+1
            ind=xt(k,:);
            ddAddt(ind(1),ind(2),i,i)=-A(ind(1),ind(2));
        end
    
        for j=(i+1):nt
            ind=xt(j,:);
            ddAddt(ind(1),ind(2),i,j)=-A(ind(1),ind(2))*tant(j)/tant(i);
            for k=(j+1):(nt+1)
                ind=xt(k,:);
                ddAddt(ind(1),ind(2),i,j)=A(ind(1),ind(2))/tant(i)/tant(j);
            end
            ddAddt(:,:,j,i)=ddAddt(:,:,i,j);
        end
    end
    
    %
    dAdf=zeros(d,d,nf); % derivative of A wrt phi: dA/df(m)
    dAdfdf=zeros(d,d,nf,nf); % second derivative of A wrt phi
    for m=1:nf
        ind=ft(m,:);
        dAdf(ind(1),ind(2),m)=-1i*A(ind(1),ind(2));
        dAdfdf(ind(1),ind(2),m,m)=-A(ind(1),ind(2));
    end
    
    pomt = zeros(d^2,d^2);
    for i = 1:d^2
        pomtemp = pom(:,:,i).';
        pomt(i,:) = pomtemp(:).';
    end
    %
    dpdm=zeros(num,num); % matrix of dp_ij/dm + dq/dm
    for i=1:num
        for j=1:nt     % wrt theta
    %         dpdm(i,j)=2*real(trace(cA*dAdt(:,:,j)*pom(:,:,i)));
            
            rtemp = cA*dAdt(:,:,j);
            rtemp = rtemp(:);
            
            dpdm(i,j)=2*real(pomt(i,:)*rtemp);
    
        end
        for j=nt+1:num     % wrt phi
    
    %         dpdm(i,j)=2*real(trace(cA*dAdf(:,:,j-nt)*pom(:,:,i)));
            
            rtemp = cA*dAdf(:,:,j-nt);
            rtemp = rtemp(:);
            dpdm(i,j)=2*real(pomt(i,:)*rtemp);
    
        end
    end
    
    % Jacobian determinant
    JacDet=abs(det(dpdm));
    
    %
    M=zeros(num,num,num);
    for i=1:nt     % wrt theta
        cdAdt=dAdt(:,:,i)';
        for j=1:num
            for k=1:nt
    %             M(j,k,i)=2*real(trace((cdAdt*dAdt(:,:,k)+cA*ddAddt(:,:,i,k))*pom(:,:,j)));
                
                rtemp = (cdAdt*dAdt(:,:,k)+cA*ddAddt(:,:,i,k));
                rtemp = rtemp(:);
                M(j,k,i)=2*real(pomt(j,:)*rtemp);
            end
            for k=nt+1:num
    %             M(j,k,i)=2*real(trace((cdAdt*dAdf(:,:,k-nt)+cA*dAdtdf(:,:,i,k-nt))*pom(:,:,j)));
                
                rtemp = (cdAdt*dAdf(:,:,k-nt)+cA*dAdtdf(:,:,i,k-nt));
                rtemp = rtemp(:);
                M(j,k,i)=2*real(pomt(j,:)*rtemp);
            end
        end
    end
    for i=nt+1:num     % wrt phi
        cdAdf=dAdf(:,:,i-nt)';
        for j=1:num
            for k=1:nt
    %             M(j,k,i)=2*real(trace((cdAdf*dAdt(:,:,k)+cA*dAdtdf(:,:,k,i-nt))*pom(:,:,j)));
                
                rtemp = (cdAdf*dAdt(:,:,k)+cA*dAdtdf(:,:,k,i-nt));
                rtemp = rtemp(:);
                M(j,k,i)=2*real(pomt(j,:)*rtemp);
            end
            for k=nt+1:num
    %             M(j,k,i)=2*real(trace((cdAdf*dAdf(:,:,k-nt)+cA*dAdfdf(:,:,k-nt,k-nt))*pom(:,:,j)));
                
                rtemp = (cdAdf*dAdf(:,:,k-nt)+cA*dAdfdf(:,:,k-nt,k-nt));
                rtemp = rtemp(:);
                M(j,k,i)=2*real(pomt(j,:)*rtemp);
                
            end
        end
    end
    
    %-------------------------------------------------------------------------%
    %--------------- distribution ---------------%
    p_full=zeros(1,16);
    dpdm_full=zeros(16,num);
    
    AAtemp = cA*A;
    AAtemp = AAtemp(:);
    for i=1:15
        
        
    %     p_full(i)=real(trace(cA*A*pom(:,:,i)));
        
        p_full(i)=real(pomt(i,:)*AAtemp);
        dpdm_full(i,:)=dpdm(i,:);
        
       
            
            
    end
    p_full(16)=1-sum(p_full(1:15));
    dpdm_full(16,:)=-sum(dpdm_full(1:15,:));
    
    % prior=abs(prod(p_full)); % conjugate prior, for instance
    % p_full corresponds to Q, the overall POVM
    
    %-------------------------------------------------------------------------%
    %--------------- experimental data ---------------%
    % for prior, set D=zeros(1,K) with K being the number of POVM outcomes
    
    % put in the target count
    % D=randi([0 20],1,d^2); % for demonstration, let's randomly generate a set of data 
    % D = [6; 8; 4; 5; 6; 6; 5; 5; 5; 4; 5; 4; 5; 4; 5; 4]';
    D = [10; 4; 6; 4; 7; 6; 5; 6; 5; 6; 10; 6; 5; 6; 8; 6]';
    
    
    %--------------- the likelihood ---------------%
    likelihood=1;
    for i=1:16
        likelihood=likelihood*power(p_full(i),D(i));
    end
    
    % add prior and likelihood to the Jacobian determinant
    % JacDet=JacDet*prior*likelihood;
    JacDet=JacDet*likelihood;
    
    
    %--------------- the gradient ---------------%
    % u_prior=zeros(1,num); % gradient invoked from prior distribution
    % for i=1:num 
    %     for j=1:d^2
    %         u_prior(i)=u_prior(i)+1/p_full(j)*dpdm_full(j,i); % for conjugate prior
    %         % dpdm_full is the first derivative of p_full w.r.t. the angle parameters
    %     end
    % end
    % u_prior=-1/2*u_prior; % for Jeffreys prior
    
    u_likeli=zeros(1,num); % gradient invoked from the likelihood
    for i=1:num
        for j=1:16
            u_likeli(i)=u_likeli(i)+D(j)/p_full(j)*dpdm_full(j,i);
        end
    end
    
    u=zeros(1,num); % overall gradient
    inv_dpdm=eye(num,num)/dpdm; % inv_dpdm=inv(dpdm);
    
    Mt = zeros(num,num^2);
    for i = 1:num
        Mtemp = M(:,:,i).';
        Mt(i,:) = Mtemp(:).';
    end
    
    inv_dpdmc = inv_dpdm(:);
    for i=1:num
    %     u(i)=trace(inv_dpdm*M(:,:,i))+u_prior(i)+u_likeli(i);
    %     u(i)=trace(inv_dpdm*M(:,:,i))+u_likeli(i);
        u(i)=Mt(i,:)*inv_dpdmc+u_likeli(i);
    
    
    end
   