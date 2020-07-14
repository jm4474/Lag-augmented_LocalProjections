classdef LALP
% -------------------------------------------------------------------------
% This Matlab class provides generic functions to implement bootstrap with 
% lag-augmented local projection for practical use. The functions provide
% users with a wild bootstrap with equal-tailed percentile-t confidence
% interval (Default horizons of interest: 1, 6, 12, 36, 60).
%
% To see our suggested sample of bootstrap implementation of 
% lag-augmented local projection, please run the following line:
% [results] = LALP.lalp(4, [0.5 0.9], 240, 0.5, 0.3, 0.1)
%
% This version: July 13th, 2020
% J. L. Montiel Olea & M. Plagborg-Moller
% -------------------------------------------------------------------------
    
    properties
    end
    
    methods(Static) 
    %The following functions can be directly accessed
    %using LALP.FunctionName().
        
    %% 1) Main Function
          function [results] = lalp(p, rhos, T, a, tau, alpha)
        
              %------------------------------------------------------------
              %This main function generates bootstrap point estimations and
              %bootstrap confidence intervals based on VAR(p) model.
              %------
              %INPUT:
              %a)              p: lag length                               (1 x 1) 
              %b)           rhos: values of parameter rho to loop over     (1 x 2)
              %c)              T: sample sizes to loop over                (1 x 1)
              %d)              a: parameter a                              (1 x 1)
              %e)            tau: covariance of innovations                (1 x 1)
              %f)          alpha: confidence level                         (1 x 1)
              %------
              %OUTPUT:
              %a)        results.estims: estimated impulse responses       (3-D double)
              %b)           results.ses: s.e. for impulse responses        (3-D double)
              %c)     results.cis_lower: lower bounds for percentile-t
              %                          confidence interval               (3-D double)
              %d)     results.cis_upper: upper bounds for percentile-t
              %                          confidence interval               (3-D double)
              %e)    results.cover_inds: indicator functions for coverage
              %                          probabilities                     (3-D double)
              %f) results.coverage_prob: coverage probabilities of
              %                          suggested bootstrap confidence
              %                          intervals                         (2-D double)
              %g)       results.lengths: lengths of suggested bootstrap
              %                          confidence intervals              (3-D double)
              %h) results.median_length: the median of lengths             (2-D double)
              %------------------------------------------------------------
        
              %% Settings
        
              delete(gcp('nocreate')); % close all current pools

              dgp = struct;
              
              dgp.p = p; % Lag length
              
              dgp.rhos = rhos; % Values of parameter rho to loop over
              
              dgp.Ts = T; % Sample sizes T to loop over
              
              dgp.a = a; % Parameter a
              
              dgp.tau = tau; % Parameter tau

              %% Monte Carlo simulation settings

              sim = struct;

              sim.numrep ...
                  = 100;                                % No. of repetitions

              sim.rng_seed ...
                  = 202006262;                           % Random number seed

              sim.num_workers ...
                  = 2;                                  % No. of parallel workers 
                                                        % (=0: run serial)
                                          
              % Reporting
              results_filename ...
                  = sprintf('%s%d', 'sim_var_p', dgp.p);  % File name for storing results                                                                

              %% Regression settings
              
              settings = struct;

              settings.p ...
                       = dgp.p;                        % Lag length used for estimation 
                                                   % (excluding augmented lags)

              settings.horzs ...
                       = [1 6 12 36 60];           % Horizons of interest
     
              settings.resp_var = 2;              % Index of response variable of interest

              settings.innov = 1;                 % Index of innovation of interest

              settings.no_const ...
                       = false;                    % true: omit intercept 
     
              settings.se_homosk = false;

              settings.boot_num ...
                       = 10;                      % Number of bootstrap samples

              settings.alpha ...
                       = alpha;                      % Significance level                               

              %% Suggested Specification

              specs = cell(1,1); % for all specifications, check "run_sim_var.m"
              
              % LP, lag-augmented, bootstrap: VAR
              specs{1} = {'estimator', 'lp',...
                          'lag_aug', true,...
                          'bootstrap', 'var'};
        
              %% Preliminaries

              % True A(2,:) and Sigma
              dgp.n = 2;
              dgp.A_lower = [dgp.a zeros(1,dgp.n*dgp.p-1)];
              for l=1:dgp.p
                  dgp.A_lower(2*l) = -nchoosek(dgp.p,l)*(-dgp.a)^l; % Coefficients on own lags of y_{2,t}
              end
              dgp.Sigma = [1 dgp.tau; dgp.tau 1];

              rng(sim.rng_seed);                   % Set RNG seed
                                     
              aux1   = repmat(dgp.rhos,size(dgp.Ts));
              aux2   = repmat(dgp.Ts',size(dgp.rhos))';
              dgps  = [aux1;reshape(aux2,[1,size(aux1,2)])];
              clear   aux1 aux2
              
              numdgp ...
                   = size(dgps,2);                 % No. of DGPs
              
              numhorz ...
                   = length(settings.horzs);       % No. of horizons
 
              %numspec ...
                   %= length(specs);                % No. of regression specifications
 
              numrep ...
                   = sim.numrep;                   % No. of repetitions

              % Cell array of settings shared among all specifications
              spec_shared = {'resp_var', settings.resp_var, ...
                             'innov', settings.innov, ...
                             'alpha', settings.alpha, ...
                             'no_const', settings.no_const, ...
                             'se_homosk', settings.se_homosk};

              %% Run simulations
              

              irs_true = zeros(numdgp, numhorz);
              
              estims ...
                  = zeros(numdgp, numhorz, numrep);
                                                             % Initialize matrix for results
              ses = estims;

              cis_lower ...
                  = zeros(numdgp, numhorz, numrep);          

              cis_upper = cis_lower;

              if sim.num_workers > 0
    
                  poolobj = parpool(sim.num_workers); 
                                                             % Start parallel workers
    
              end

              timer = tic; % Start timer

              for i_dgp = 1:numdgp
    
                  i_rho = dgps(1,i_dgp);
    
                  i_T   = dgps(2,i_dgp);
    
                  i_rand_seeds ...
                        = randi(2^32-1,1,numrep);           % Random number seeds to be 
                                                            % supplied to parallel workers
    
                  fprintf('%2d%s%2d %s%4.2f %s%5d\n',...
                          i_dgp,...
                          '/', ...
                          numdgp,...
                          ': rho=',i_rho,...
                          'T=', i_T);             % Display the current iteration
    
	              % True VAR parameters
                  i_A = [i_rho zeros(1,dgp.n*dgp.p-1);
                         dgp.A_lower];
    
                  % True impulse responses
                  i_ir_all = LALP.var_ir(i_A,settings.horzs);
                  irs_true(i_dgp,:) = i_ir_all(settings.resp_var,settings.innov,:);
    
                  parfor(i=1:numrep, sim.num_workers) % For each repetition...
        
                      rng(i_rand_seeds(i));       % Set RNG seed
    
                      % Simulate VAR(p) data
        
                      i_Y = LALP.var_sim(i_A, zeros(dgp.n,1), mvnrnd(zeros(1,dgp.n),dgp.Sigma,i_T), zeros(dgp.p,dgp.n)); % Data series (with y_0=...=y_{1-p}=0)
        
                      i_estims ...
                          = zeros(1,numhorz);
        
                      i_ses ...
                          = i_estims;
        
                      i_cis_lower ...
                          = nan(1,numhorz);  
        
                      i_cis_upper ...
                          = i_cis_lower;
        
                      % Run all estimation procedures and compute delta method CIs
            
                          % Estimate (delta method CIs are not used in this file)
                          [i_estims(:),...
                           i_ses(:),   ...
                           i_cis_dm,   ...
                           i_cis_boot] = LALP.ir_estim(i_Y, settings.p, settings.horzs, spec_shared{:}, specs{1}{:});
            
            
                          % Bootstrap confidence intervals
                          i_cis_lower(:) = i_cis_boot(1,:);
                          i_cis_upper(:) = i_cis_boot(2,:);
                          
                      %end
                      
                      % Store all results for this repetition
                      estims(i_dgp,:,i) = i_estims;
                      ses(i_dgp,:,i) = i_ses;
                      cis_lower(i_dgp,:,i) = i_cis_lower(:,:);
                      cis_upper(i_dgp,:,i) = i_cis_upper(:,:);
                      
                      if mod(i, ceil(numrep/10)) == 0
                          fprintf('%6d%s\n', i/numrep*100, '%');
                      end
                      
                  end
    
              end

              sim.elapsed_time = toc(timer);

              disp('Elapsed time (min):');
              disp(sim.elapsed_time/60);

              if sim.num_workers > 0
                  delete(poolobj); % Stop parallel workers
              end


              %% Compute coverage and median length

              dgp.dgps ...
                  = dgps;

              dgp.irs_true = irs_true;

              % Store results
              results ...
                  = struct;

              results.estims ...
                  = estims; 
              
              results.ses ...
                  = ses;
              
              results.cis_lower ...
                  = cis_lower;
              
              results.cis_upper ...
                  = cis_upper;
              
              % Coverage
              irs_true_reshape ...
                  = dgp.irs_true;
              
              results.cover_inds ...
                  = (irs_true_reshape >= results.cis_lower ...
                    & irs_true_reshape <= results.cis_upper) + 0; 
                    % Coverage indicator (1 or 0)
      
              results.cover_inds(isnan(results.cis_lower)) ...
                  = nan; 
                    % Set indicator to missing if no CI is recorded
      
              results.coverage_prob ...
                  = mean(results.cover_inds, 3); 
                    % Coverage probability

              % Length
              results.lengths ...
                  = results.cis_upper-results.cis_lower; % Lengths
              
              results.median_length ...
                  = median(results.lengths, 3); % Median length
          end
    
    %% 2) Local Projection
          function [ir, ir_varcov, betahat, betahat_varcov, res, X]...
                    = lp(Y,num_lags,horz,resp_ind,se_setting,no_const)

              %------------------------------------------------------------
              %This function generates local projection coefficients and 
              %variance-covariance matrices with lagged controls.
              %------
              %INPUT:
              %a)             Y: data matrix                               (T x n)
              %b)     n um_lags: number of lags of Y to control for
              %                 (in addition to the contemporaneous
              %                 regressor)                                 (1 x 1)
              %c)          horz: horizon of interest                       (1 x 1)
              %d)      resp_ind: indices of response variables
              %                 Y_t of interest                            (1 x m)
              %e)    se_setting: EITHER bool: if true, homoskedastic s.e.;
              %                               if false, EHW s.e.
              %                  OR function handle: function that returns
              %                  HAC/HAR sandwich matrix
              %f)      no_const: true: omit intercept                      (bool)
              %------
              %OUTPUT:
              %a)             ir: estimated impulse responses
              %                   at select horizons                       (m x n)
              %b)      ir_varcov: var-cov of vec(ir)                       (mn x mn)
              %c)        betahat: full vector of estimated regression
              %                   coefficients                             ((np+n+1) x 1)
              %d) betahat_varcov: var-cov of betahat                       ((np+n+1) x (np+n+1)) 
              %e)            res: residuals                                ((T-p) x m)
              %f)              X: covariate data matrix including intercept((T-p) x (np+n+1))
              %------------------------------------------------------------
    
              % Covariate matrix
              Y_lag = lagmatrix(Y,0:num_lags);
              X = Y_lag(num_lags+1:end-horz,:);
    
              % Local projection
              [betahat, betahat_varcov, res, X] = LALP.linreg(Y(num_lags+horz+1:end,resp_ind),X,se_setting,no_const);
    
              n = size(Y,2);
              m = length(resp_ind);
              ir = betahat(:,1:n);
              ir_varcov = betahat_varcov(1:m*n,1:m*n);

          end
    
    %% 3) System Linear Regression
          function [betahat, varcov, res, X_expand] = linreg(Y, X, se_setting, no_const)

              %------------------------------------------------------------
              %This function computes system linear regression
              %Y_t = beta*X_t + epsilon_t
              %where dim(Y_t) = n, dim(beta) = n x k.
              %------
              %INPUT:
              %a)          Y: dependent variable data matrix               (T x n)
              %b)          X: covariate data matrix                        (T x k)
              %c) se_setting: EITHER bool: if true, homoskedastic s.e.;
              %                            if false, EHW s.e.
              %               OR function handle: 
              %               function that returns HAC/HAR sandwich matrix
              %d)   no_const: true: omit intercept                         (bool)
              %------
              %OUTPUT:
              %a)    betahat: estimated coefficients                       (n x (k+1))
              %b)     varcov: var-cov matrix of vec(betahat)               (n(k+1) x n(k+1))
              %c)        res: residual matrix                              (T x n)
              %d)   X_expand: expanded covariate data matrix with intercept(T x (k+1))
              %------------------------------------------------------------
              
              [T,n] = size(Y);
    
              % Include intercept if desired
              X_expand = [X ones(T,1-no_const)];
              k = size(X_expand,2);
    
              % OLS
              betahat = (X_expand\Y)';
    
              % Standard errors
    
              if nargout > 1
        
                  res = Y-X_expand*betahat';
                  XpX = X_expand'*X_expand;
                  scores = kron(X_expand,ones(1,n)).*repmat(res,1,k);

                  if islogical(se_setting)
                      if se_setting % If homoskedastic s.e.
                          varcov = kron(inv(XpX),(res'*res)/T);
                      else % If EHW s.e.
                          varcov = LALP.kron_fast(inv(XpX),LALP.kron_fast(inv(XpX),scores'*scores,0)',0)'; % EHW var-cov matrix
                      end
                  else % If HAC/HAR s.e.
                      varcov = LALP.kron_fast(inv(XpX), LALP.kron_fast(inv(XpX),se_setting(scores),0)',0)'; % HAC/HAR var-cov matrix
                  end

                  varcov = T/(T-k)*varcov; % Finite sample adjustment as in Stata
    
              end

          end
    
    %% 4) Estimation of VAR(p) model
    
          function [irs, irs_varcov, Ahat_estim, res_estim]...
                    = var_ir_estim(Y, p, p_estim, horzs, bias_corr, homosk, no_const)
    
              %------------------------------------------------------------
              %This function computes VAR(p) least-squares estimates
              %and delta method s.e., allowing for lag augmentation.
              %------
              %INPUT:
              %a)         Y: data vector  (T x n)
              %b)         p: lag length used 
              %              for impulse response computations             (1 x 1)       
              %c)   p_estim: lag length used for estimation (p_estim >= p) (1 x 1)
              %d)     horzs: horizons of interest                          (H x 1)
              %e) bias_corr: true: apply analytical bias correction        (bool)
              %f)    homosk: true: homoskedastic s.e., false: EHW s.e.     (bool)
              %g)  no_const: true: omit intercept                          (bool)
              %------
              %OUTPUT:
              %a)       irs: estimated impulse responses
              %              Theta_h at select horizons                    (n x n x H)
              %b)irs_varcov: var-cov matrices of vec(Theta_h)
              %              at select horizons                            (n^2 x n^2 x H)
              %c)Ahat_estim: VAR coefficient estimates [A_1,...,A_p]
              %              (possibly bias-corrected,
              %              possibly including intercept as last column)  (n x np)
              %d) res_estim: estimation residuals                          ((T-p_estim) x n)
              %------------------------------------------------------------
    
              [T,n] = size(Y);
    
              % One-step forecasting regression of Y_{t+1} on (Y_t, ..., Y_{t-p_estim+1})
              [~,~,Ahat_estim,Ahat_estim_varcov,res_estim] = LALP.lp(Y,p_estim-1,1,1:n,homosk,no_const);
    
              % If bias correction is desired...
              if bias_corr
                  Sigmahat = (res_estim'*res_estim)/(size(res_estim,1)-n*p_estim-1+no_const); % Residual variance estimate
                  Ahat_estim(:,1:end-1+no_const) = LALP.var_biascorr(Ahat_estim(:,1:end-1+no_const), Sigmahat, T);
              end
    
              % Only use first p VAR coefficient matrices to compute impulse responses
              Ahat = Ahat_estim(:,1:n*p);
              Ahat_varcov = Ahat_estim_varcov(1:n^2*p,1:n^2*p);
    
              if nargout==1
                  irs = LALP.var_ir(Ahat,horzs); % Compute impulse responses
              else
                  [irs, jacob] = LALP.var_ir(Ahat,horzs); % Compute impulse responses and Jacobian
                  nh = length(horzs);
                  irs_varcov = zeros(n^2,n^2,nh);
                  for h=1:nh
                      irs_varcov(:,:,h) = jacob(:,:,h)*Ahat_varcov*jacob(:,:,h)'; % Var-cov for impulse response matrix at this horizon
                  end
              end

          end
    
    
    %% 5) Bias-adjusted VAR coefficients
    
          function A_corr = var_biascorr(A, Sigmahat, T)

              %------------------------------------------------------------
              %This function adjusts bias for VAR(p) estimator
              %based on Pope (JTSA 1990), equation 9.
              %------
              %INPUT:
              %a)        A: original VAR(p) coefficient estimates
              %             [A_1,...,A_p]                                  (n x np)
              %b) Sigmahat: estimate of VAR(p) innovation variance         (n x n)
              %c)        T: sample size                                    (1 x 1)
              %------
              %OUTPUT:
              %a)   A_corr: bias-corrected coefficient estimates           (n x np)
              %------------------------------------------------------------
    
              %Set up companion form: X_t = A*X_{t-1} + Z_t, where dim(X_t)=n*p
              [n,np] = size(A);
              A_comp = [A; eye(np-n), zeros(np-n,n)];
    
              if max(abs(eig(A_comp)))>1 
                 %If original point estimate is outside stationary region, do not bias correct
                 A_corr = A;
              return;
              end
    
              G = blkdiag(Sigmahat, zeros(np-n)); % Var(Z_t)
              Gamma0 = reshape((eye(np^2)-kron(A_comp,A_comp))\G(:),np,np); % Var(X_t)
    
              %Bias correction formula
              aux = inv(eye(np)-A_comp')+(A_comp')/(eye(np)-A_comp'*A_comp');
              lambdas = eig(A_comp);
              for lamb = lambdas'
                  aux = aux + lamb*inv(eye(np)-lamb*A_comp');
              end
              b = G*(aux/Gamma0); % Scaled negative bias
              A_corr = A_comp + b/T; % Bias-corrected companion form coefficients
    
              %If corrected estimate is outside stationary region, 
              %reduce bias correction little by little
              %(as recommended by Kilian & Lütkepohl, 2017, ch. 12)
              delta = 1;
              while max(abs(eig(A_corr)))>1 && delta>0
                  delta = delta - 0.01;
                  A_corr = A_comp + delta*b/T;
              end
    
              %Return bias-corrected VAR(p) coefficients
              A_corr = A_corr(1:n,:);

          end
          
    %% 6) VAR(p) IRF and Jacobian
    
          function [irs, jacob] = var_ir(A,horzs)

              %--------------------------------------------------------------
              %This function computes VAR(p) reduced-form impulse responses
              %and Jacobian based on given VAR coefficients.
              %------
              %INPUT:
              %a)     A: VAR coefficient matrices (A_1, ..., A_p)          (n x np)
              %b)  horz: horizons of interest                              (1 x H)
              %------
              %OUTPUT:
              %a)   irs: reduced-form impulse responses Theta_h
              %          at select horizons                                (n x n x H)
              %b) jacob: Jacobian of vec(Theta_h) at select
              %          horizons wrt. vec(A)                              (n^2 x (n^2*p) x H)
              %------------------------------------------------------------
    
    
              % Dimensions
              nh = length(horzs);
              maxh = max(horzs);
              [n,np] = size(A);
              p = np/n;
    
              irs = zeros(n,n,nh); % Will contain impulse responses at select horizons
              jacob = zeros(n^2,n^2*p,nh); % Will contain Jacobian of impulse responses at select horizons wrt. vec(A)
    
              ir_p = [eye(n); zeros(n*(p-1),n)]; % Will contain last p impulse responses, stacked vertically
              jacob_p = zeros(n^2,n^2*p,p); % Will contain last p values of the Jacobian of vec(Theta_h) wrt. vec(A)
    
              for h=1:maxh % Loop through horizons
        
                  the_A = A(:,1:n*min(h,p));
                  the_past_ir = ir_p(1:n*min(h,p),:);
                  the_ir = the_A*the_past_ir; % Impulse response at horizon h
                  ir_p = [the_ir; ir_p(1:end-n,:)]; % Shift forward in time
        
                  the_ind = find(horzs==h,1);
                  if ~isempty(the_ind)
                      irs(:,:,the_ind) = the_ir; % Store impulse response at select horizons
                  end
        
                  if nargout>1
            
                      % Jacobian at horizon h via chain rule
                      the_jacob_p = zeros(n^2,n^2*p);
                      the_jacob_p(:,1:n^2*min(h,p)) = kron(the_past_ir',eye(n));
                      for l=1:min(h,p)
                          the_jacob_p(:,1:n^2*min(h,p)) = the_jacob_p(:,1:n^2*min(h,p)) + LALP.kron_fast(A(:,(l-1)*n+1:l*n),jacob_p(:,1:n^2*min(h,p),l),1);
                      end
            
                      % Shift forward in time
                      jacob_p(:,:,2:end) = jacob_p(:,:,1:end-1);
                      jacob_p(:,:,1) = the_jacob_p;
            
                      if ~isempty(the_ind)
                          jacob(:,:,the_ind) = the_jacob_p; % Store Jacobian at select horizons
                      end
            
                  end
        
              end

          end
    
    %% 7) VAR(p) Simulation
          function Y = var_sim(A, c, U, Y_init)

              %------------------------------------------------------------
              %This function simulates VAR(p) data
              %given innovations and initial conditions.
              %------
              %INPUT:
              %a)       A: VAR coefficients [A_1,...,A_p]                  (n x np)
              %b)       c: Intercepts                                      (n x 1)
              %c)       U: Innovations data matrix                         (T x n?
              %d)  Y_init: Initial conditions [y_1,...,y_p]'
              %            (only last p rows will be used)                ?p x n?
              %------
              %OUTPUT:
              %a)       Y: Simulated data matrix [y_1,...,y_T]',
              %            discarding initial obs.                         (T x n)
              %------------------------------------------------------------
              
              % Dimensions
              [T,n] = size(U);
              p = size(A,2)/n;
              
              % Iterate
              Y = zeros(T+p,n);
              Y(1:p,:) = Y_init(end-p+1:end,:); % Initial conditions
              for t=p+1:T+p
                  Y(t,:) = c' + reshape(Y(t-1:-1:t-p,:)',1,n*p)*A' + U(t-p,:);
              end
              Y = Y(p+1:end,:); % Discard initial conditions

          end
    
    %% 8) Wild Bootstrap based on VAR(p) Model
    
          function Y_boot = var_boot(A, res, Y, p, wild, no_const)

              %------------------------------------------------------------
              %This function generates wild bootstrap sample by iterating
              %on the VAR(p) model.
              %------
              %INPUT:
              %a)       A: VAR(p) coefficient matrices [A_1,...,A_p]       (n x np)
              %            OPTIONAL: length could exceed n*p, in which case 
              %                      last column equals intercept
              %                NOTE: only A(:,1:n*p) and A(:,end) will be used
              %b)     res: residuals                                       (T_res x n)
              %c)       Y: data vector                                     (T x n)
              %d)       p: lag length                                      (1 x 1)
              %e)    wild: true: wild bootstrap,
              %            false: homoskedastic bootstrap                  (bool)
              %f)no_const: true: exclude intercept from bootstrap samples  (bool)
              %------
              %OUTPUT:
              %a)  Y_boot: bootstrap sample                                (T x n)
              %------------------------------------------------------------

              % Dimensions
              [T,n] = size(Y);
              T_res = size(res,1); 
              % Effective sample size for residuals (need not equal T-p)
    
              % Draw block of initial T-T_res obs. from real data
              ind_init = randi(T_res+1);
              Y_init = Y(ind_init:ind_init+T-T_res-1,:);
    
              % Intercept (if applicable)
              c = zeros(n,1);
              if ~no_const
                  c = A(:,end);
              end

              % Draw residuals
              if wild % Wild bootstrap (default)
                  res_boot = randn(T_res,1).*res; 
              else    % I.i.d. bootstrap
                  res_boot = res(randi(T_res,T_res,1),:);
              end
    
              % Generate VAR(p) data, with residuals and initial conditions as above
              Y_boot = [Y_init; LALP.var_sim(A(:,1:n*p), c, res_boot, Y_init)];

          end
          
    %% 9) Bootstrap Confidence Intervals
      
          function cis = boot_ci(pseudo_truth, estim_data, ses_data, estims_boot, ses_boot, alpha)
          
              %----------------------------------------------------------------
              %This function generates bootstrap confidence intervals.
              %(The default setting generates Hall percentile-t confidence
              % intervals).
              %------
              %INPUT:
              %a) pseudo_truth: pseudo-true parameter vector
              %                 used to generate bootstrap samples             (1 x k)
              %b)   estim_data: parameter estimates in real data               (1 x k)
              %c)     ses_data: s.e. in real data                              (1 x k)
              %d)  estims_boot: bootstrapped parameter estimates               (B x k)
              %e)     ses_boot: bootstrapped s.e.                              (B x k)
              %f)        alpha: significance level                             (1 x 1)
              %------
              %OUTPUT:
              %a)          cis: bootstrap intervals
              %                 (1st index: lower and upper limit,
              %                  2nd index: parameter,
              %                  3rd index: type of interval - 
              %                             Efron, Hall, or Hall percentile-t) (2 x k x 3)
              %----------------------------------------------------------------
        
        
              % Bootstrap quantiles
              estim_boot_quants ...
                  = quantile(estims_boot, [alpha/2 1-alpha/2], 1); 
                    % Bootstrap quantiles of estimates
              
          tstats_boot ...
              = (estims_boot - pseudo_truth)./ses_boot; 
                % Bootstrap t-stats
          
          tstat_boot_quants ...
              = quantile(tstats_boot, [1-alpha/2 alpha/2], 1); 
                % Bootstrap quantiles of t-stats
    
          % Compute confidence intervals
          ci_efron ...
              = estim_boot_quants; 
              % Efron
        
          ci_hall ...
              = estim_data + pseudo_truth - estim_boot_quants([2 1],:); 
              % Hall
        
          ci_hall_t ...
              = estim_data - ses_data.*tstat_boot_quants; 
              % Hall percentile-t
    
          % Reshape
          cis = reshape([ci_hall_t], 2, []); 
          % (If all types of confidence intervals are needed, please
          %  uncomment the next line)
          %cis = reshape([ci_efron ci_hall ci_hall_t], 2, [], 3);
          

          end
          
    %% 10) Returning Impulse Responses (VAR & LP)
          function [irs, ses] = var_select(irs_all, irs_all_varcov, resp_var, nu)

              %------------------------------------------------------------
              %This function computes impulse responses of interest
              %along with s.e..
              %------------------------------------------------------------

              irs = nu'*permute(irs_all(resp_var,:,:), [2 3 1]);
              
              % Standard errors
              if nargout>1
                  [n,~,nh] = size(irs_all);
                  the_eye = eye(n);
                  aux = kron(nu',the_eye(resp_var,:));
                  ses = zeros(1,nh);
                  for h=1:nh
                      ses(h) = sqrt(aux*irs_all_varcov(:,:,h)*aux');
                  end
              end

          end
          
          function [ir, se] = lp_select(irs_all, irs_all_varcov, nu)

              % Local projection: Return impulse response of interest along with s.e.
              ir = irs_all*nu;
              se = sqrt(nu'*irs_all_varcov*nu);

          end
    
    %% 11) Fast Kronecker Product Function
          
          function C = kron_fast(A,B,type)
              
              %------------------------------------------------------------
              %This function generates efficient computation of
              %C = kron(A, eye(n))*B (type=0)
              %or
              %C = kron(eye(n), A)*B (type=1)
              %------
              %From "MATLAB array manipulation tips and tricks"
              %by Peter J. Acklam (2002), section 10.8
              %------------------------------------------------------------
    
              % Matrix dimensions
              [p,q] = size(A);
              [qn,m] = size(B);
              n = qn/q;
    
              % Compute
              switch type
                  case 0
                      C = reshape(reshape(B.', [n*m q])*A.', [m p*n]).';
                  case 1
                      C = reshape(A*reshape(B, [q n*m]), [p*n m]);
              end

          end    
          
    %% 12) Impulse Response Estimations
          function [irs, ses, cis_dm, cis_boot] = ir_estim(Y, p, horzs, varargin)

              %------------------------------------------------------------
              %This is the wrapper function for AR or LP estimation of
              %impulse responses delta method and bootstrap confidence intervals
              %------
              %INPUT:
              %See "Parse inputs" section
              %------
              %OUTPUT:
              %a)      irs: estimated impulse responses at select horizons (1 x H)
              %b)      ses: s.e. for impulse responses                     (1 x H)
              %c)   cis_dm: lower and upper limits of
              %             delta method confidence intervals              (2 x H)
              %d) cis_boot: lower and upper limits of bootstrap
              %             confidence intervals
              %             (3rd index: type of interval, either Efron,
              %             Hall, or Hall percentile-t)                    (2 x H x 3)
              %------------------------------------------------------------
    
              %% Parse inputs
    
              ip = inputParser;
    
              % Required inputs
              addRequired(ip, 'Y', @isnumeric);
                  % T x n     data vector
              addRequired(ip, 'p', @isnumeric);
                  % 1 x 1     lag length (not counting any augmentation)
              addRequired(ip, 'horzs', @isnumeric);
                  % 1 x H     impulse response horizons of interest
    
              % Optional inputs
              addParameter(ip, 'resp_var', 1, @isnumeric);
                  % Index of response variable of interest (default: first variable)
              addParameter(ip, 'innov', 1, @isnumeric);
                  % Index of innovation of interest, or n x 1 vector with linear combination of innovations (default: first innovation)
              addParameter(ip, 'estimator', 'lp', @ischar);
                  % Estimator type, either 'var' or 'lp' (default: local projection)
              addParameter(ip, 'alpha', 0.05, @isnumeric);
                  % Significance level (default: 0.05)
              addParameter(ip, 'lag_aug', false, @islogical);
                  % Lag-augment? (default: no)
              addParameter(ip, 'bias_corr', true, @islogical);
                  % Bias-correct VAR estimates? (default: yes)
              addParameter(ip, 'se_homosk', false, @islogical);
                  % Homoskedastic standard errors/bootstrap? (default: no)
              addParameter(ip, 'no_const', false, @islogical);
                  % Omit intercept? (default: no)
              addParameter(ip, 'har', [], @(x) isa(x, 'function_handle') || isempty(x));
                  % HAR estimator as function of data and bandwidth, empty if no HAR (default: no HAR)
              addParameter(ip, 'har_bw', [], @(x) isa(x, 'function_handle') || isempty(x));
                  % HAR bandwidth as function of sample size
              addParameter(ip, 'har_cv', [], @(x) isa(x, 'function_handle') || isempty(x));
                  % HAR critical value as function of bandwidth, empty if normal critical value
              addParameter(ip, 'bootstrap', [], @(x) ischar(x) || isempty(x));
                  % Bootstrap type, either 'var', 'resid', or 'pair'; empty if delta method (default: delta method)
              addParameter(ip, 'boot_num', 1000, @isnumeric);
                  % Bootstrap iterations (default: 1000)
              addParameter(ip, 'boot_lag_aug', false, @islogical);
                  % Lag-augment in bootstrap DGP? Only relevant for 'var' bootstrap (default: no)
    
              parse(ip, Y, p, horzs, varargin{:});
    
    
              %% Preliminaries
    
              [T,n] = size(Y);              % Dimensions
    
              nh ...
                = length(horzs); % Number of horizons
              
              cvs ...
                = repmat(norminv(1-ip.Results.alpha/2),...
                        1, nh);      % Default critical values: normal
          
              cis_boot ...
                = nan(2,nh,3);       % Initializes NaN array
  
              % Determine linear combination of innovations
              if length(ip.Results.innov)==1
                  the_eye = eye(n);
                  nu = the_eye(:,ip.Results.innov); % Unit vector
              else
                  nu = ip.Results.innov(:); % User-specified vector
              end
   
              %% Point estimates and var-cov
    
              if strcmp(ip.Results.estimator, 'var') % VAR
        
                  % VAR impulse responses
                  [irs_all, irs_all_varcov] = LALP.var_ir_estim(Y, ...
                                                           p,...
                                                           p+ip.Results.lag_aug,...
                                                           horzs, ...
                                                           ip.Results.bias_corr,...
                                                           ip.Results.se_homosk,...
                                                           ip.Results.no_const);
        
                  % Impulse responses of interest and s.e.
                  [irs, ses] = LALP.var_select(irs_all, irs_all_varcov, ip.Results.resp_var, nu);
        
              elseif strcmp(ip.Results.estimator, 'lp') % LP
        
                  irs = zeros(1,nh);
        
                  ses = zeros(1,nh);
                  
                  betahat ...
                      = cell(1,nh);
        
                  res = cell(1,nh);
        
                  X   = cell(1,nh);
        
                  for h=1:nh % For each horizon...
            
                      the_horz = horzs(h); % Horizon
            
                      % Determine s.e. setting and c.v.
            
                      if ~isempty(ip.Results.har) % HAR
                          the_bw ...
                              = ip.Results.har_bw(T-p-ip.Results.lag_aug-the_horz); 
                              % Bandwidth, determined by effective sample size
                    
                          the_se_setting ...
                              = @(X) ip.Results.har(X,the_bw); 
                              % HAR function
                 
                          if ~isempty(ip.Results.har_cv)
                              cvs(h) ...
                                  = ip.Results.har_cv(the_bw); % Critical value
                          end
                
                      else % EHW/homoskedastic
                
                          the_se_setting = ip.Results.se_homosk; % Indicator for whether homosk. or EHW
            
                      end
            
                      % LP regression
                      [the_irs_all, the_irs_all_varcov, betahat{h},...
                      ~, res{h}, X{h}] = LALP.lp(Y, ...
                                            p-1+ip.Results.lag_aug,...                                                                                          
                                            the_horz,...      
                                            ip.Results.resp_var,...
                                            the_se_setting,...
                                            ip.Results.no_const);
                      [irs(h),ses(h)] = LALP.lp_select(the_irs_all, the_irs_all_varcov, nu);
            
                  end
                  
              end
    
              % If only point estimates and standard errors are requested, stop
              if nargout<=2
                  return;
              end
    
    
              %% Delta method confidence intervals
              
              cis_dm = irs + [-1; 1]*(cvs.*ses);
    
    
              %% Bootstrap confidence intervals
              
              if ~isempty(ip.Results.bootstrap)
        
                  estims_boot = zeros(ip.Results.boot_num,nh);
                  ses_boot = estims_boot;
        
                  if strcmp(ip.Results.bootstrap, 'var') % Recursive VAR bootstrap specifications
            
                      % VAR coefficient estimates that define bootstrap DGP
            
                      [~, ~, Ahat_var, res_var] = LALP.var_ir_estim(Y, ...
                                                            p,...
                                                            p+ip.Results.boot_lag_aug,...
                                                            horzs, ...
                                                            ip.Results.bias_corr,...
                                                            ip.Results.se_homosk,...
                                                            ip.Results.no_const);
            
                      pseudo_truth = LALP.var_select(LALP.var_ir(Ahat_var(:,1:n*(p+ip.Results.boot_lag_aug)), horzs), [], ip.Results.resp_var, nu); % Pseudo-true impulse responses in bootstrap DGP
            
                      for b=1:ip.Results.boot_num

                          % Generate bootstrap sample based on (possibly lag-augmented) AR estimates
                          Y_boot = LALP.var_boot(Ahat_var, res_var, Y, p+ip.Results.boot_lag_aug, ip.Results.se_homosk, ip.Results.no_const);

                          % Estimate on bootstrap sample
                          [estims_boot(b,:),ses_boot(b,:)] = LALP.ir_estim(Y_boot, p, horzs, varargin{:});

                      end

                  else % Linear regression bootstrap specifications
            
                      pseudo_truth = irs;
                      
                      for h=1:nh % Treat each horizon separately

                          for b=1:ip.Results.boot_num

                              % Generate bootstrap sample
                              [Y_boot, X_boot] = linreg_boot(betahat{h}',res{h},X{h},ip.Results.bootstrap,ip.Results.se_homosk);

                              % Run OLS on bootstrap sample
                              [the_linreg_betahat, the_linreg_varcov] = LALP.linreg(Y_boot,X_boot,ip.Results.se_homosk,true); % Don't add extra intercept
                              [estims_boot(b,h),ses_boot(b,h)] = LALP.lp_select(the_linreg_betahat(1:n),the_linreg_varcov(1:n,1:n),nu);

                          end

                      end
        
                  end
        
                  % Compute bootstrap confidence intervals
                  cis_boot = LALP.boot_ci(pseudo_truth, irs, ses, estims_boot, ses_boot, ip.Results.alpha);
        
              end

          end
    end

end