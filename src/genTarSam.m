function [prob_points_accepted,ar] = genTarSam(N_total,pop,kappa_w,n_w,lambp,kappa_s,n_s,lambso,lambs)
    % 
    % This generates the target sample acoording to the target poterior with a specified reference sample
    % 
    % 
    % Input
    % --------------------------------------------------------------------------
    % N_total : int
    %   # of sample points generated
    % pop : array of integers
    %   the power of probabilities for the posterior
    % kappa_w : real
    %   percentage of Wishart reference sample in the total reference sample
    % n_w : int
    %   # of columns of the Psi matrix used in construction of the Wishart sample
    % lambp : real
    %   the 'distance' between the origin and the Wishart ref peak
    % kappa_s : real
    %   percentage of shifted reference sample in the total reference sample
    % n_s : int
    %   # of columns of the Psi matrix used in construction of the shifted sample
    % lambso : real
    %   the 'distance' between the origin and the ref peak before linear shift
    % lambs : real
    %   to which place the ref sample is shifted, 1 being MLE, 0 being origin
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % prob_points_accepted : 2d array of complex double (dimension : m^2 x N)
    %   desired sample points in the probability space, distributed according to the target posterior 
    % ar : real
    %   acceptance rate
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see main.m
    % 
    % 
    % 
    % 
    
    if N_total > 5e7
        N_total_ori = N_total;
        N_trial = N_total / 1e7;
        N_total = 1e7;
        save_dir = '..\temp\';
    else
        N_trial = 1;
    end

    n_qubit = log2(length(pop)) / 2;
    kappa = 1 - kappa_w - kappa_s;
    
    if not(floor(n_qubit) == n_qubit)
        error("the length of the posterior is not a power of an integer number of qubits");
    end
    
    pom = buildNTetraPOM(n_qubit); % the default measurement is multiple tetrahedron
    m = 2 ^ n_qubit;
    
    % split the large sample into smaller runs with 1e5 points each
    if N_total > 1e5 
        N_run = ceil(N_total/1e5);
        N_ref_single = 1e5;
    else
        N_run = 1;
        N_ref_single = N_total;
    end
    
    %--- MLE
    ref_freq = pop / sum(pop); % relative frequency
    % linear inversion from probability to rho
    rho_peak = prob2Rho(ref_freq,pom);     
    % to check if the linearly inverted rho is physical
    v = eig(rho_peak);

    if sum(v<0)
        fprintf('linear inversion not physical, mle running \n');
        rhomle = qse_apg(pom,ref_freq');
        fprintf('mle done \n');
        rho_peak = lambp*rhomle+(1-lambp)*eye(m)/m;   
        ppeak = rho2Prob(rhomle,pom)';
        %-- o. peak for the shifted ref
        rho_peak_so = rhomle*lambso+(1-lambso)*eye(m)/m;
        rho_peak_s = rhomle*lambs+(1-lambs)*eye(m)/m;
    else
        rho_peak_so = rho_peak*lambso+(1-lambso)*eye(m)/m;
        rho_peak_s = rho_peak*lambs+(1-lambs)*eye(m)/m; 
        ppeak = rho2Prob(rho_peak,pom)';
    end  
    
    % log of likelihood of the maximum likelihood estimator
    logLpeak = sum(pop.*log(ppeak));

    N_w = ceil(N_ref_single*kappa_w);
    N_s = ceil(N_ref_single*kappa_s);
    N_uni = N_ref_single - N_w - N_s;

    log_ratio_max = zeros(1,N_trial);
    
    for t_dx = 1 : N_trial
        log_ratio_runs = zeros(N_ref_single,N_run);
        prob_points_runs = zeros(m^2,N_ref_single,N_run);

        parfor nrunk = 1 : N_run
            %-- obtaining the ref sam
            [rhos,~,Sigma_w,Sigma_s] = genRefSam(n_w,n_s,N_w,N_s,N_ref_single,m,rho_peak,rho_peak_so,rho_peak_s);
            %prob_points/rhos returns the sample in the p/rho-space
            %Sigma are the convariance matrices
            %indunph records the index of the points that could only be from 
            %the Wishart but not the linear shift part so that the resampling
            %ratio can be calculated correctly


            %-- calculation of the reference probability of the points
            ref_prob = calRefProb(rhos,kappa,kappa_w,n_w,Sigma_w,n_s,Sigma_s,rho_peak_so,rho_peak_s);
            
            %-- calculation of the target probability of the points
            [log_tar_prob,prob_points] = calTarProb(rhos,pop,pom);

            %-- the acceptance ratio in log is "lal"
            log_ref_prob = log(ref_prob);
            log_ratio = log_tar_prob - log_ref_prob;

            log_ratio_runs(:,nrunk) = log_ratio;
            prob_points_runs(:,:,nrunk) = prob_points;
        end

        log_ratio_trial = log_ratio_runs(:);
        [log_ratio_max(t_dx),~] = max(log_ratio_trial);
        prob_points_trial = reshape(prob_points_runs,[m^2,N_total]);

        if N_trial > 1
            % If the original N_total are splitted into different 
            save([save_dir 'temp_data',num2str(t_dx) '.mat'],'prob_points_trial','log_ratio_trial','-v7.3');
        else
            ldraw = log(rand(size(log_ratio_trial)));
            acc = ldraw < (log_ratio_trial-log_ratio_max);
            prob_points_accepted = prob_points_trial(:,acc);
            ar = size(prob_points_accepted,2) / N_total;
        end
    end

    if N_trial > 1
        [log_ratio_max,~] = max(log_ratio_max);
        prob_points_accepted = [];
        for t_dx = 1 : N_trial
            load([save_dir 'temp_data',num2str(t_dx) '.mat']);
            ldraw = log(rand(size(log_ratio_trial)));
            acc = ldraw < (log_ratio_trial-log_ratio_max);
            prob_points_accepted = [prob_points_accepted,prob_points_trial(:,acc)];
            delete([save_dir 'temp_data',num2str(t_dx) '.mat']);
        end
        ar = size(prob_points_accepted,2) / N_total_ori;
    end


