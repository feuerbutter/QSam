function ref_prob = calRefProb(rhos,kappa,kappa_w,n_w,Sigma,varargin)
    % 
    % This calculates the reference probability according to Eqn.(16)
    % 
    % 
    % Input
    % --------------------------------------------------------------------------
    % rhos : 3d array of comnplex double
    %   sample points in the state space
    % kappa : double
    %   percentage of uniform sample
    % kappa_w : double
    %   percentage of Wishart sample
    % n_w : int   
    %   # of columns of the Psi matrix used in construction of the Wishart sample
    % Sigma : array of double
    %   the covariance matrix required to produce the peaked Wishart sample
    % varargin : 4 additional optional parameters to include shifted reference sample
    %   n_s = varargin{1};
    %   Sigma_s = varargin{2};
    %   rho_peak_shift_o = varargin{3};
    %   rho_peak_shift = varargin{4};
    % 
    % 
    % Output
    % --------------------------------------------------------------------------
    % ref_prob : 1d array of real
    %   the reference probability of the sample points
    % 
    % How to call
    % --------------------------------------------------------------------------
    % see genTarSam.m
    % 
    % 
    % 

    kappa_s = 1 - kappa - kappa_w;

    N = size(rhos,3);
    m = size(rhos,2);
    ref_prob = zeros(1,N);

    inv_Sigma = inv(Sigma);
    norm_wis = calWishartNorm(n_w,m,inv_Sigma);
    norm_uni = calWishartNorm(m,m,eye(m));

    if isempty(varargin)
        for r_dx = 1 : N 
            ref_prob(r_dx) = kappa_w * exp(calWishartProb(rhos(:,:,r_dx),n_w,m,inv_Sigma)+norm_wis-norm_uni) + kappa;
        end
    else
        n_s = varargin{1};
        Sigma_s = varargin{2};
        rho_peak_shift_o = varargin{3};
        rho_peak_shift = varargin{4};
        inv_Sigma_s = inv(Sigma_s);
        norm_shf = calWishartNorm(n_s,m,inv_Sigma_s);
        delrho = rho_peak_shift - rho_peak_shift_o;
        for r_dx = 1 : N 
            rho = rhos(:,:,r_dx);
            if max(rho(:)) == 0
                ref_prob(r_dx) = 0;
            else
                if  min(real(eig(rho-delrho))) < 0
                    % if (rho-delrho) is not physical, that means it cannot be from the shifted sample, then only the Wishart and uniform probability needs to be calculated, while in the else condition, all three probabilities need to be computed
                    ref_prob(r_dx) = kappa_w * exp(calWishartProb(rho,n_w,m,inv_Sigma)+norm_wis-norm_uni) + kappa;
                else
                    ref_prob(r_dx) = kappa_w * exp(calWishartProb(rho,n_w,m,inv_Sigma)+norm_wis-norm_uni) + ...
                    kappa_s * exp(calWishartProb(rho-delrho,n_s,m,inv_Sigma_s)+norm_shf-norm_uni) + ...
                    kappa;
                end
            end
        end
    end
    