clear;
addpath('functions/');

% Monte Carlo study of VAR(p) inference procedures


%% Settings

dgp = struct;

dgp.rhos = [0.5 0.9];

dgp.Ts = [240 2400];


%% Monte Carlo simulation settings

sim = struct;

sim.numrep ...
    = 1e3;                                % No. of repetitions

sim.rng_seed ...
    = 20200625;                           % Random number seed

sim.num_workers ...
    = 4;                                  % No. of parallel workers 
                                          % (=0: run serial)


%% Regression settings

settings = struct;

settings.p ...
         = 1;                        % Lag length used for estimation 
                                     % (excluding augmented lags)

settings.horzs ...
         = [1 6 12 36 60];           % Horizons of interest
     
settings.resp_var = 2;
settings.innov = 1;

settings.no_const ...
         = false;                    % true: omit intercept 
     
settings.se_homosk = false;

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

specs = cell(3,1);           % Specifications for the simulations

specs{1} = {'estimator', 'ar',...
            'lag_aug', false,....
            'bias_corr', false};

specs{2} = {'estimator', 'lp',...
            'lag_aug', true};

specs{3} = {'estimator', 'lp',...
            'lag_aug', true,...
            'har', settings.har};

        
%% Preliminaries

rng(sim.rng_seed);                   % Set RNG seed

%dgps = combvec(dgp.rhos, dgp.Ts);   % DGPs 
                                     % combvec Requires Matlab 2019
                                     % The following line works in previous
                                     % versions
                                     % Consider changing the name of dgps
                                     % to persistence_horizon_combination
                                     % to avoid confusion with the
                                     % structure dgp.
                                     
aux1   = repmat(dgp.rhos,size(dgp.Ts));
aux2   = repmat(dgp.Ts',size(dgp.rhos))';
dgps  = [aux1;reshape(aux2,[1,size(aux1,2)])];
clear   aux1 aux2

numdgp ...
     = size(dgps,2);                 % No. of DGPs

numhorz ...
     = length(settings.horzs);       % No. of horizons
 
numspec ...
     = length(specs);                % No. of regression specifications
 
numrep ...
     = sim.numrep;                   % No. of repetitions

% Cell array of settings shared among all specifications
spec_shared = {'resp_var', settings.resp_var, ...
               'innov', settings.innov, ...
               'alpha', settings.alpha, ...
               'no_const', settings.no_const, ...
               'se_homosk', settings.se_homosk, ...
               'har_bw', settings.har_bw, ...
               'har_cv', settings.har_cv};


%% Run simulations

irs_true = zeros(numdgp, numhorz);

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
    
    % True impulse responses
    i_beta = [i_rho 0; 0.5 0.5];
    i_ir_all = var_ir(i_beta,settings.horzs);
    irs_true(i_dgp,:) = i_ir_all(settings.resp_var,settings.innov,:);
    i_n = size(i_beta,1);
    
    parfor(i=1:numrep, sim.num_workers) % For each repetition...
%     for i=1:numrep
        
        rng(i_rand_seeds(i));       % Set RNG seed
    
        % Simulate VAR(p) data
        
        i_Y = var_sim(i_beta, randn(i_T,i_n), zeros(settings.p,i_n)); % Data series (with y_0=...=y_{1-p}=0)
        
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
           
