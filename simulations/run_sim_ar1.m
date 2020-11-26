clear;
addpath('../functions/');

% Monte Carlo study in AR(1) model

% Estimation procedures considered:
% 1. AR, w/ or w/o lag augmentation, w/ or w/o analytical bias correction
% 2. LP, w/ or w/o lag augmentation

% Inference procedures considered:
% a. Delta method (homoskedastic, EHW, or HAR)
% b. Recursive AR bootstrap (homoskedastic or wild, lag-augmented or not)
% c. OLS residual bootstrap (homoskedastic or wild)
% d. OLS pairs/nonparametric bootstrap

% Bootstrap confidence intervals considered:
% i. Efron percentile interval
% ii. Hall percentile interval
% iii. Hall percentile-t interval


% MPM 2019-11-22


%% Data Generating Process (DGP)
dgp = struct;

dgp.type ...
    = 'iid';                    % ='iid': i.i.d. innovations.
                                % But s.e./bootstraps allow for heterosk.; 
                                % ='arch': ARCH innovations;
                                % ='homosk': i.i.d. innovations and homosk.
                                % s.e./bootstraps;

dgp.rhos ...
    = [0 0.5 0.9 0.95 1];       % AR(1) parameters to consider

dgp.Ts ...
    = [240 480 2400];           % Sample sizes to consider

% GARCH(1,1) innovation process
% u_t = sigma_t*epsilon_t
% sigma_t^2 = alpha_0 + alpha_1*u_{t-1}^2 + beta_1*sigma_{t-1}^2
% epsilon_t ~ N(0,1)
if strcmp(dgp.type, 'arch')
    
    dgp.garch_alpha_1 = 0.7;
    
else
    
    dgp.garch_alpha_1 = 0;
    
end

dgp.garch_beta_1 = 0;

dgp.garch_alpha_0 = 1 - dgp.garch_alpha_1 - dgp.garch_beta_1; 
                                % So that E[sigma_t^2]=1
 
%% Monte Carlo simulation settings

sim = struct;

sim.numrep ...
    = 5e3;                                % No. of repetitions

sim.rng_seed ...
    = 202007251;                           % Random number seed

sim.num_workers ...
    = 4;                                  % No. of parallel workers 
                                          % (=0: run serial)

% Reporting
results_filename ...
    = sprintf('%s%s', 'sim_ar1_', dgp.type);  % File name for storing results                                                                


%% Regression settings

settings = struct;

settings.p ...
         = 1;                        % Lag length used for estimation 
                                     % (excluding augmented lags)
     
settings.horzs ...
         = [1 6 12 36 60];           % Horizons of interest
     
settings.no_const ...
         = false;                    % true: omit intercept 
     
if strcmp(dgp.type, 'homosk')
    
    settings.se_homosk = true;       % false: EHW s.e./wild bootstrap, 
                                     % true: homoskedastic s.e./bootstrap 
                                     % (doesn't apply to HAR procedure)
else
    
    settings.se_homosk = false;
    
end

settings.boot_num ...
         = 2e3;                      % Number of bootstrap samples
                                     
settings.alpha ...
         = 0.1;                      % Significance level

settings.har ....
         = @(Y,bw) ewc(Y,bw);        % HAR estimator

settings.har_bw ...
         = @(T) round(0.4*T.^(2/3)); % HAR bandwidth

settings.har_cv ...
         = @(bw) tinv(1-settings.alpha/2,bw); ...
                                     % HAR critical value


%% List of specifications for the simulations

specs = cell(9,1);           % Specifications for the simulations

specs{1} = {'estimator', 'var',... %Estimation Method: "plug-in" ar
            'lag_aug', false,....  %With or Without Lag-augmentation
            'bias_corr', false};   %With or Without Biascorrection
                                   % If 'bootstrap' option is not
                                   % specified, standard errors are
                                   % based on the Eicker-Huber-White
                                   % formula.
                                   % Heteroskedasticity and
                                   % Autocorrelation robust standard
                                   % errors can be incorporated with
                                   % the 'har' option (see spec 8)
                                        
specs{2} = {'estimator', 'var',...
            'lag_aug', false,...
            'bias_corr', true};
              
specs{3} = {'estimator', 'var',...
            'lag_aug', true,...
            'bias_corr', false};
              
specs{4} = {'estimator', 'var',...
            'lag_aug', true,...
            'bias_corr', true,...
            'bootstrap', 'var',...
            'boot_lag_aug', true};

specs{5} = {'estimator', 'var',...
            'lag_aug', true,...
            'bias_corr', true,...
            'bootstrap', 'var',...
            'boot_lag_aug', false};

specs{6} = {'estimator', 'lp', ...
            'lag_aug', false,...
            'har', settings.har,...
            'bootstrap', 'var',...
            'boot_lag_aug', false};
              
specs{7} = {'estimator', 'lp',...
            'lag_aug', true,...
            'bootstrap', 'var',...
            'boot_lag_aug', false};

specs{8} = {'estimator', 'lp',...
            'lag_aug', true,...
            'bootstrap', 'resid'};
              
specs{9} = {'estimator', 'lp',...
            'lag_aug', true, ...
            'bootstrap', 'pair'};

%% Preliminaries

rng(sim.rng_seed, 'twister');                   % Set RNG seed

% Combinations of DGP parameters
aux1   = repmat(dgp.rhos,size(dgp.Ts));
aux2   = repmat(dgp.Ts',size(dgp.rhos))';
dgps  = [aux1;reshape(aux2,[1,size(aux1,2)])];
clear   aux1 aux2;

numdgp ...
     = size(dgps,2);                 % No. of DGPs
 
numhorz ...
     = length(settings.horzs);       % No. of horizons
 
numspec ...
     = length(specs);                % No. of regression specifications
 
numrep ...
     = sim.numrep;                   % No. of repetitions

% Cell array of settings shared among all specifications
spec_shared = {'alpha', settings.alpha, ...
               'no_const', settings.no_const, ...
               'se_homosk', settings.se_homosk, ...
               'boot_num', settings.boot_num, ...
               'har_bw', settings.har_bw, ...
               'har_cv', settings.har_cv};


%% Run simulations

estims ...
    = zeros(numdgp, numspec, numhorz, numrep);
                                     % Initialize matrix for results
ses = estims;

cis_lower ...
    = zeros(numdgp, numspec, numhorz, 4, numrep); 
                                     % 4th index = type of CI: 
                                     % i) delta method, 
                                     % ii) Efron, 
                                     % iii) Hall, 
                                     % iv) Hall percentile-t
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
          = randi(2^32-1,1,numrep); % Random number seeds to be 
                                    % supplied to parallel workers
    
    fprintf('%2d%s%2d %s%4.2f %s%5d\n',...
            i_dgp,...
            '/', ...
            numdgp,...
            ': rho=',i_rho,...
            'T=', i_T);             % Display the current iteration
    
    parfor(i=1:numrep, sim.num_workers) % For each repetition...
%     for i=1:numrep
        
        rng(i_rand_seeds(i), 'twister');       % Set RNG seed
    
        % Simulate AR(1) data with GARCH innovations
        i_U = sim_garch(dgp.garch_alpha_0,...
                        dgp.garch_alpha_1,...
                        dgp.garch_beta_1,...
                        i_T);
        
        i_Y = filter(1, [1 -i_rho], i_U); % Data series (with y_0=0)
        
        i_estims ...
            = zeros(numspec, numhorz);
        
        i_ses ...
            = i_estims;
        
        i_cis_lower ...
            = nan(numspec, numhorz, 4);  % 4 refers to the 4 types of CIs
        
        i_cis_upper ...
            = i_cis_lower;
        
        % Run all estimation procedures and compute delta method CIs
        for j=1:numspec
            
            % Estimate
            [i_estims(j,:),...
             i_ses(j,:),   ...
             i_cis_dm,   ...
             i_cis_boot] = ir_estim(i_Y, settings.p, settings.horzs, spec_shared{:}, specs{j}{:});
            
            % Delta method confidence interval
            i_cis_lower(j,:,1) = i_cis_dm(1,:);
            i_cis_upper(j,:,1) = i_cis_dm(2,:);
            
            % Bootstrap confidence intervals
            i_cis_lower(j,:,2:end) = i_cis_boot(1,:,:);
            i_cis_upper(j,:,2:end) = i_cis_boot(2,:,:);
            
        end
        
        % Store all results for this repetition
        estims(i_dgp,:,:,i) = i_estims;
        ses(i_dgp,:,:,i) = i_ses;
        cis_lower(i_dgp,:,:,:,i) = i_cis_lower;
        cis_upper(i_dgp,:,:,:,i) = i_cis_upper;
        
        if mod(i, ceil(numrep/10)) == 0
            fprintf('%6d%s\n', round(i/numrep*100), '%');
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

dgp.irs_true ...
    = dgps(1,:)'.^settings.horzs; % True impulse responses

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
    = permute(dgp.irs_true, [1 3 2]);

results.cover_inds ...
    = (irs_true_reshape >= results.cis_lower ...
      & irs_true_reshape <= results.cis_upper) + 0; 
      % Coverage indicator (1 or 0)
      
results.cover_inds(isnan(results.cis_lower)) ...
    = nan; 
      % Set indicator to missing if no CI is recorded
      
results.coverage_prob ...
    = mean(results.cover_inds, 5); 
      % Coverage probability

% Length
results.lengths ...
    = results.cis_upper-results.cis_lower; % Lengths

results.median_length ...
    = median(results.lengths, 5); % Median length

%% Save results
status = mkdir('results');
save(strcat('results/',results_filename, '.mat'), 'dgp', 'specs', 'settings', 'sim', 'results');
