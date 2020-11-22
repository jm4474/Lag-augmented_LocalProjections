clear;
addpath('../functions/');

% Monte Carlo study of VAR inference procedures
% calibrated to Gertler & Karadi (AEJ:Macro 2015) data

% MPM 2020-11-20


%% Data settings

data = struct;

% Data files from Gertler & Karadi (2015)
data.file_var = '../examples/gk_data/VAR_data.csv';      % VAR series
data.file_iv = '../examples/gk_data/factor_data.csv';    % External instrument series

% Variables
data.vars = {'logcpi', 'logip', 'ff', 'ebp', 'ff4_tc'};  % Variables in VAR


%% Monte Carlo simulation settings

sim = struct;

sim.numrep ...
    = 1e3;                                % No. of repetitions

sim.rng_seed ...
    = 202011211;                           % Random number seed

sim.num_workers ...
    = 4;                                  % No. of parallel workers 
                                          % (=0: run serial)
                                          
% Reporting
results_filename ...
    = sprintf('%s%d', 'sim_var_calib');  % File name for storing results                                                                


%% Regression settings

settings = struct;

settings.p = 12;                        % Lag length used for estimation 
                                     % (excluding augmented lags)

settings.horzs ...
         = 1:48;           % Horizons of interest
     
settings.resp_vars = [1 2 4]; % Indices of response variables of interest

settings.innov = 5; % Index of innovation of interest

settings.no_const ...
         = false;                    % true: omit intercept 
     
settings.se_homosk = false;

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

specs = cell(4,1);           % Specifications for the simulations

% VAR, non-augmented
specs{1} = {'estimator', 'var',...
            'lag_aug', false,...
            'bootstrap', 'var'};

% VAR, lag-augmented
specs{2} = {'estimator', 'var',...
            'lag_aug', true,...
            'bootstrap', 'var', ...
            'boot_lag_aug', true};

% LP, non-augmented, bootstrap: VAR
specs{3} = {'estimator', 'lp',...
            'lag_aug', false,...
            'har', settings.har,...
            'bootstrap', 'var'};

% LP, lag-augmented, bootstrap: VAR
specs{4} = {'estimator', 'lp',...
            'lag_aug', true,...
            'bootstrap', 'var'};

        
%% Estimate VAR in data

% Load data
data.dat = innerjoin(readtable(data.file_var), readtable(data.file_iv));
data.Y = data.dat{:,data.vars}; % Select variables
data.Y = data.Y(~any(isnan(data.Y),2),:); % Remove missing
data.Y = detrend(data.Y,0); % Remove mean
[data.T,data.n] = size(data.Y); % Sample size

% Estimate VAR and true IRFs
[data.irs, ~, data.A, data.res] ...
    = var_ir_estim(data.Y, settings.p, settings.p, settings.horzs, ...
                   false, false, false);
               
% Residual var-cov matrix
data.Sigma = (data.res'*data.res)/(size(data.res,1)-size(data.A,2));

% Display largest eigenvalues
data.A_comp ...
    = [data.A(:,1:end-1);
      eye(data.n*(settings.p-1)) zeros(data.n*(settings.p-1),data.n)];
                                                    % Companion matrix
disp('Five largest VAR eigenvalues (absolute values)');
disp(abs(eigs(data.A_comp,5))');


%% Preliminaries

% True parameters
dgp = struct;
dgp.A = data.A(:,1:end-1);
dgp.Sigma = data.Sigma;
dgp.irs = data.irs;
dgp.n = data.n;
dgp.T = data.T;

rng(sim.rng_seed, 'twister');                   % Set RNG seed

resp_vars = settings.resp_vars;     % Indices of response variables

numvar ...
     = length(settings.resp_vars);                 % No. of response variables

numhorz ...
     = length(settings.horzs);       % No. of horizons
 
numspec ...
     = length(specs);                % No. of regression specifications
 
numrep ...
     = sim.numrep;                   % No. of repetitions

% Cell array of settings shared among all specifications
spec_shared = {'p', settings.p, ...
               'innov', settings.innov, ...
               'alpha', settings.alpha, ...
               'no_const', settings.no_const, ...
               'se_homosk', settings.se_homosk, ...
               'boot_num', settings.boot_num, ...
               'har_bw', settings.har_bw, ...
               'har_cv', settings.har_cv};


%% Run simulations

estims ...
    = zeros(numvar, numspec, numhorz, numrep);
                                     % Initialize matrix for results
ses = estims;

cis_lower ...
    = zeros(numvar, numspec, numhorz, 4, numrep); 
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

rand_seeds ...
      = randi(2^32-1,1,numrep); % Random number seeds to be 
                                % supplied to parallel workers

timer = tic; % Start timer

parfor(i=1:numrep, sim.num_workers) % For each repetition...
% for i=1:numrep
    
    rng(rand_seeds(i), 'twister');       % Set RNG seed

    % Simulate VAR(p) data

    i_Y = var_sim(dgp.A, zeros(dgp.n,1), mvnrnd(zeros(1,dgp.n),dgp.Sigma,dgp.T), zeros(settings.p,dgp.n)); % Data series (with y_0=...=y_{1-p}=0)

    i_estims ...
        = zeros(numvar, numspec, numhorz);

    i_ses ...
        = i_estims;

    i_cis_lower ...
        = nan(numvar, numspec, numhorz, 4);  % 4 refers to the 4 types of CIs

    i_cis_upper ...
        = i_cis_lower;

    for i_var = 1:numvar % For each response variable of interest...

        i_resp_var = resp_vars(i_var);

        % Run all estimation procedures and compute delta method CIs
        for j=1:numspec

            % Estimate
            [i_estims(i_var,j,:),...
             i_ses(i_var,j,:),   ...
             i_cis_dm,   ...
             i_cis_boot] = ir_estim(i_Y, settings.p, settings.horzs, ...
                                    'resp_var', i_resp_var, ...
                                    spec_shared{:}, specs{j}{:});

            % Delta method confidence interval
            i_cis_lower(i_var,j,:,1) = i_cis_dm(1,:);
            i_cis_upper(i_var,j,:,1) = i_cis_dm(2,:);

            % Bootstrap confidence intervals
            i_cis_lower(i_var,j,:,2:end) = i_cis_boot(1,:,:);
            i_cis_upper(i_var,j,:,2:end) = i_cis_boot(2,:,:);

        end

    end

    % Store all results for this repetition
    estims(:,:,:,i) = i_estims;
    ses(:,:,:,i) = i_ses;
    cis_lower(:,:,:,:,i) = i_cis_lower;
    cis_upper(:,:,:,:,i) = i_cis_upper;

    if mod(i, ceil(numrep/50)) == 0
        fprintf('%6d%s\n', round(i/numrep*100), '%');
    end
    
end

sim.elapsed_time = toc(timer);

disp('Elapsed time (min):');
disp(sim.elapsed_time/60);

if sim.num_workers > 0
    delete(poolobj); % Stop parallel workers
end


%% Compute coverage and median length

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
irs_true ...
    = dgp.irs(settings.resp_vars,settings.innov,:);

results.cover_inds ...
    = (irs_true >= results.cis_lower ...
      & irs_true <= results.cis_upper) + 0; 
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
save(strcat('results/',results_filename, '.mat'), 'data', 'dgp', 'specs', 'settings', 'sim', 'results');