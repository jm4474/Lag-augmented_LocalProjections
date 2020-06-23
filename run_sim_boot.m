clear;
addpath('functions/');

% Monte Carlo study of AR(1) inference procedures

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


%% Settings

% DGP
dgp = struct;
dgp.type = 'iid';               % ='iid': i.i.d. innovations but s.e./bootstraps allow for heterosk.; ='arch': ARCH innovations; ='homosk': i.i.d. innovations and homosk. s.e./bootstraps;
dgp.rhos = [0 0.5 0.9 0.95 1];  % AR(1) parameters to consider
dgp.Ts = [240 480 2400];        % Sample sizes to consider

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
dgp.garch_alpha_0 = 1 - dgp.garch_alpha_1 - dgp.garch_beta_1; % Normalize E[sigma_t^2]=1

% List of inference specifications
boot_ar_spec_nonaug = struct('estimator', 'ar', 'lag_aug', false, 'biascorr', true); % Non-augmented specification used for generating bootstrap samples
boot_ar_spec_lagaug = struct('estimator', 'ar', 'lag_aug', true,  'biascorr', true); % Lag-augmented specification used for generating bootstrap samples
specs = cell(9,1);
specs{1} = struct('estimator', 'ar', 'lag_aug', false, 'biascorr', false);
specs{2} = struct('estimator', 'ar', 'lag_aug', false, 'biascorr', true);
specs{3} = struct('estimator', 'ar', 'lag_aug', true,  'biascorr', false);
specs{4} = struct('estimator', 'ar', 'lag_aug', true,  'biascorr', true, 'bootstrap', 'ar', 'bootstrap_spec', boot_ar_spec_lagaug);
specs{5} = struct('estimator', 'ar', 'lag_aug', true,  'biascorr', true, 'bootstrap', 'ar', 'bootstrap_spec', boot_ar_spec_nonaug);
specs{6} = struct('estimator', 'lp', 'lag_aug', false, 'har', true,      'bootstrap', 'ar', 'bootstrap_spec', boot_ar_spec_nonaug);
specs{7} = struct('estimator', 'lp', 'lag_aug', true,  'har', false,     'bootstrap', 'ar', 'bootstrap_spec', boot_ar_spec_nonaug);
specs{8} = struct('estimator', 'lp', 'lag_aug', true,  'har', false,     'bootstrap', 'resid');
specs{9} = struct('estimator', 'lp', 'lag_aug', true,  'har', false,     'bootstrap', 'pair');

% Regression settings
settings = struct;
settings.p = 1;                                         % Lag length used for estimation (excluding augmented lags)
settings.horzs = [1 6 12 36 60];                        % Horizons of interest
settings.noconst = false;                               % true: omit intercept in all regressions
if strcmp(dgp.type, 'homosk')
    settings.se_homosk = true;                          % false: EHW s.e./wild bootstrap, true: homoskedastic s.e./bootstrap (doesn't apply to HAR procedure)
else
    settings.se_homosk = false;
end
settings.numboot = 2e3;                                 % Number of bootstrap samples
settings.boot_ar_restrict = []; %@(x) min(max(x,-1),1); % Function that restricts coefficients to stationary region when generating AR bootstrap samples (if empty [], don't restrict)
settings.alpha = 0.1;                                   % Significance level
settings.har_bw = @(T) round(0.4*T.^(2/3));             % HAR bandwidth
settings.har_fct = @(Y,bw) ewc(Y,bw);                   % HAR estimator
settings.har_cv = @(bw) tinv(1-settings.alpha/2,bw);    % HAR critical value

% Monte Carlo simulation settings
sim = struct;
sim.numrep = 5e3;               % No. of repetitions
sim.rng_seed = 20191203;        % Random number seed
sim.num_workers = 4;            % No. of parallel workers (=0: run serial)

% Reporting
results_filename = sprintf('%s%s', 'sim_', dgp.type);  % File name for storing results


%% Preliminaries

rng(sim.rng_seed); % Set RNG seed

dgps = combvec(dgp.rhos, dgp.Ts); % DGPs
numdgp = size(dgps,2); % No. of DGPs
numhorz = length(settings.horzs); % No. of horizons
numspec = length(specs); % No. of regression specifications
numrep = sim.numrep; % No. of repetitions


%% Run simulations

estims = zeros(numdgp, numspec, numhorz, numrep);
ses = estims;
cis_lower = zeros(numdgp, numspec, numhorz, 4, numrep); % 4th index = type of CI: i) delta method, ii) Efron, iii) Hall, iv) Hall percentile-t
cis_upper = cis_lower;

if sim.num_workers > 0
    poolobj = parpool(sim.num_workers); % Start parallel workers
end

timer = tic; % Start timer

for s=1:numdgp
    
    the_rho = dgps(1,s);
    the_T = dgps(2,s);
    the_rand_seeds = randi(2^32-1,1,numrep); % Random number seeds to be supplied to parallel workers
    
    fprintf('%2d%s%2d %s%4.2f %s%5d\n', s, '/', numdgp, ': rho=', the_rho, 'T=', the_T);
    
    parfor(i=1:numrep, sim.num_workers) % For each repetition...
%     for i=1:numrep
        
        rng(the_rand_seeds(i)); % Set RNG seed
    
        % Simulate AR(1) data with GARCH innovations
        the_U = sim_garch(dgp.garch_alpha_0, dgp.garch_alpha_1, dgp.garch_beta_1, the_T); % Innovations
        the_Y = filter(1, [1 -the_rho], the_U); % Data series (with y_0=0)
        
        the_estims = zeros(numspec, numhorz);
        the_ses = the_estims;
        the_cis_lower = nan(numspec, numhorz, 4);
        the_cis_upper = the_cis_lower;
        
        % Run all estimation procedures and compute delta method CIs
        for j=1:numspec
            
            % Estimate
            [the_estims(j,:),the_ses(j,:),the_cis_dm,the_cis_boot] = ir_estim(the_Y, specs{j}, settings);
            
            % Delta method confidence interval
            the_cis_lower(j,:,1) = the_cis_dm(1,:);
            the_cis_upper(j,:,1) = the_cis_dm(2,:);
            
            % Bootstrap confidence intervals
            the_cis_lower(j,:,2:end) = the_cis_boot(1,:,:);
            the_cis_upper(j,:,2:end) = the_cis_boot(2,:,:);
            
        end
        
        % Store all results for this repetition
        estims(s,:,:,i) = the_estims;
        ses(s,:,:,i) = the_ses;
        cis_lower(s,:,:,:,i) = the_cis_lower;
        cis_upper(s,:,:,:,i) = the_cis_upper;
        
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

dgp.dgps = dgps;
dgp.irs_true = dgps(1,:)'.^settings.horzs; % True impulse responses

% Store results
results = struct;
results.estims = estims;
results.ses = ses;
results.cis_lower = cis_lower;
results.cis_upper = cis_upper;

% Coverage
irs_true_reshape = permute(dgp.irs_true, [1 3 2]);
results.cover_inds = (irs_true_reshape >= results.cis_lower & irs_true_reshape <= results.cis_upper) + 0; % Coverage indicator (1 or 0)
results.cover_inds(isnan(results.cis_lower)) = nan; % Set indicator to missing if no CI is recorded
results.coverage_prob = mean(results.cover_inds, 5); % Coverage probability

% Length
results.lengths = results.cis_upper-results.cis_lower; % Lengths
results.median_length = median(results.lengths, 5); % Median length


%% Save results

save(strcat(results_filename, '.mat'), 'dgp', 'specs', 'settings', 'sim', 'results');


