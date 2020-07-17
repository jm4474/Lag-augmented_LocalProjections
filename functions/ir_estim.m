function [irs, ses, cis_dm, cis_boot] = ir_estim(Y, p, horzs, varargin)

    % Wrapper function for AR or LP estimation of impulse responses
    % Delta method and bootstrap confidence intervals
    
    % Inputs: see below
    
    % Outputs:
    % irs       1 x H       estimated impulse responses at select horizons
    % ses       1 x H       s.e. for impulse responses
    % cis_dm    2 x H       lower and upper limits of delta method confidence intervals
    % cis_boot  2 x H x 3   lower and upper limits of bootstrap confidence intervals (3rd index: type of interval, either Efron, Hall, or Hall percentile-t)
    
    
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
    addParameter(ip, 'boot_workers', 0, @isnumeric);
        % Number of parallel workers used for bootstrapping (default: 0, meaning no parallel computation)
    addParameter(ip, 'verbose', false, @islogical);
        % Print progress to screen when bootstrapping? (default: no)
    
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
        [irs_all, irs_all_varcov] = var_ir_estim(Y, ...
                                                 p,...
                                                 p+ip.Results.lag_aug,...
                                                 horzs, ...
                                                 ip.Results.bias_corr,...
                                                 ip.Results.se_homosk,...
                                                 ip.Results.no_const);
        
        % Impulse responses of interest and s.e.
        [irs, ses] = var_select(irs_all, irs_all_varcov, ip.Results.resp_var, nu);
        
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
            ~, res{h}, X{h}] = lp(Y, ...
                                  p-1+ip.Results.lag_aug,...                                                                                          
                                  the_horz,...      
                                  ip.Results.resp_var,...
                                  the_se_setting,...
                                  ip.Results.no_const);
            [irs(h),ses(h)] = lp_select(the_irs_all, the_irs_all_varcov, nu);
            
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
            
            [irs_var, ~, Ahat_var, res_var] = var_ir_estim(Y, ...
                                                  p+ip.Results.boot_lag_aug,... % Compute impulse responses for full lag-augmented model, if applicable
                                                  p+ip.Results.boot_lag_aug,...
                                                  horzs, ...
                                                  ip.Results.bias_corr,...
                                                  ip.Results.se_homosk,...
                                                  ip.Results.no_const);
            
            pseudo_truth = var_select(irs_var, [], ip.Results.resp_var, nu); % Pseudo-true impulse responses in bootstrap DGP
            
%             for b=1:ip.Results.boot_num
            parfor(b=1:ip.Results.boot_num, ip.Results.boot_workers)

                % Generate bootstrap sample based on (possibly lag-augmented) VAR estimates
                Y_boot = var_boot(Ahat_var, res_var, Y, p+ip.Results.boot_lag_aug, ip.Results.se_homosk, ip.Results.no_const);

                % Estimate on bootstrap sample
                [estims_boot(b,:),ses_boot(b,:)] = ir_estim(Y_boot, p, horzs, varargin{:});
                
                % Print progress
                print_prog(b, ip.Results.boot_num, ip.Results.verbose);

            end

        else % Linear regression bootstrap specifications
            
            pseudo_truth = irs;
            
%             for b=1:ip.Results.boot_num
            parfor(b=1:ip.Results.boot_num, ip.Results.boot_workers)
  
                for h=1:nh % Treat each horizon separately
                    
                    % Generate bootstrap sample
                    [Y_boot, X_boot] = linreg_boot(betahat{h}',res{h},X{h},ip.Results.bootstrap,ip.Results.se_homosk);

                    % Run OLS on bootstrap sample
                    [the_linreg_betahat, the_linreg_varcov] = linreg(Y_boot,X_boot,ip.Results.se_homosk,true); % Don't add extra intercept
                    [estims_boot(b,h),ses_boot(b,h)] = lp_select(the_linreg_betahat(1:n),the_linreg_varcov(1:n,1:n),nu);

                end
                
                % Print progress
                print_prog(b, ip.Results.boot_num, ip.Results.verbose);

            end
        
        end
        
        % Compute bootstrap confidence intervals
        cis_boot = boot_ci(pseudo_truth, irs, ses, estims_boot, ses_boot, ip.Results.alpha);
        
    end

end

function [irs, ses] = var_select(irs_all, irs_all_varcov, resp_var, nu)

    % VAR: Return impulse responses of interest along with s.e.

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

function print_prog(b, n, verbose)
   
    % Print progress to screen periodicially
    
    if ~verbose
        return;
    end
    
    if mod(b,ceil(n/10))==0
        fprintf('%3d%s\n', 100*b/n, '%');
    end
    
end