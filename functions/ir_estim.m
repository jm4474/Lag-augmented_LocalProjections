function [irs, ses, cis_dm, cis_boot, betahat, res] = ir_estim(Y, spec, settings)

    % Wrapper function for AR or LP estimation of impulse responses
    % Delta method and bootstrap confidence intervals
    
    % Inputs:
    % Y         T x 1   data vector
    % spec      struct  estimator specification
    % settings  struct  estimation settings
    % varargin          OPTIONAL: estimator specification used to generate AR bootstrap samples
    
    % Outputs:
    % irs       H x 1       estimated impulse responses
    % ses       H x 1       s.e. for impulse responses
    % cis_dm    2 x H       lower and upper limits of delta method confidence intervals
    % cis_boot  2 x H x 3   lower and upper limits of bootstrap confidence intervals (3rd index: type of interval, either Efron, Hall, or Hall percentile-t)
    % betahat               coefficient estimates in regression(s)
    % res                   residuals in regression(s)
    
    
    %% Preliminaries
    
    T = length(Y);              % Sample size
    
    numhorz ...
      = length(settings.horzs); % Number of horizons
    
    cvs ...
      = repmat(norminv(1-settings.alpha/2),...
              1, numhorz);      % Default critical values: normal
          
    cis_boot ...
      = nan(2,numhorz,3);       % Initializes NaN array.
    
    
    %% Point estimates and s.e.
    
    if strcmp(spec.estimator, 'ar') % AR
        
        % AR impulse responses
        [irs, ses, ~, betahat, res] = ar_ir_estim(Y, ...
                                                  settings.p,...
                                                  settings.p+spec.lag_aug,...
                                                  settings.horzs, ...
                                                  spec.biascorr,...
                                                  settings.se_homosk,...
                                                  settings.noconst);
        
    elseif strcmp(spec.estimator, 'lp') % LP
        
        irs = zeros(numhorz,1);
        
        ses = zeros(numhorz,1);
        
        betahat ...
            = cell(numhorz,1);
        
        res = cell(numhorz,1);
        
        X   = cell(numhorz,1);
        
        for h=1:numhorz % For each horizon...
            
            the_horz = settings.horzs(h); % Horizon
            
            % Determine s.e. setting and c.v.
            
            if spec.har % HAR
                the_bw ...
                    = settings.har_bw(T-settings.p-spec.lag_aug-the_horz); 
                    % Bandwidth, determined by effective sample size
                    
                the_se_setting ...
                    = @(X) settings.har_fct(X,the_bw); 
                    % HAR function
                    
                cvs(h) ...
                    = settings.har_cv(the_bw); % Critical value
                
            else % EHW/homoskedastic
                
                the_se_setting = settings.se_homosk; % Indicator for whether homosk. or EHW
            
            end
            
            % LP regression
            %settings p is the dimension of the AR model
            %spec.lag_aug is whether or not the LP is lag augmented 
            [irs(h), ses(h), betahat{h},...
            ~, res{h}, X{h}] = lp(Y, ...
                                  settings.p-1+spec.lag_aug,...                                                                                          
                                  the_horz,...                  
                                  the_se_setting,...
                                  settings.noconst);
            
        end
        
    end
    
    
    %% Delta method confidence intervals
    
    cis_dm = irs' + [-1; 1]*(cvs.*ses');
    
    
    %% Bootstrap confidence intervals
    
    if isfield(spec, 'bootstrap') && ~isempty(spec.bootstrap)
        
        estims_boot = zeros(settings.numboot,numhorz);
        ses_boot = estims_boot;
        
        if strcmp(spec.bootstrap, 'ar') % Recursive AR bootstrap specifications
            
            spec_boot_ar = spec.bootstrap_spec; 
                         % Specification for generating bootstrap samples
                         
            spec_boot_ar.bootstrap ...
                         = []; % Don't do another bootstrap when estimating bootstrap DGP
            
            % AR coefficient estimates that define bootstrap DGP
            
            %[~, ~, ~, ~, betahat_ar, res_ar] ...
            %            = ir_estim(Y, spec_boot_ar, settings);
            [~, ~, ~, betahat_ar, res_ar] = ar_ir_estim(Y, ...
                                                  settings.p,...
                                                  settings.p+spec_boot_ar.lag_aug,...
                                                  settings.horzs, ...
                                                  spec_boot_ar.biascorr,...
                                                  settings.se_homosk,...
                                                  settings.noconst);
                                              
            if isfield(settings, 'boot_ar_restrict') && ~isempty(settings.boot_ar_restrict)
                
                betahat_ar(1:settings.p) = settings.boot_ar_restrict(betahat_ar(1:settings.p)); % Restrict coefficient estimates to stationary region
            
            end
            
            pseudo_truth = ar_ir(betahat_ar(1:settings.p), settings.horzs); % Pseudo-true impulse responses in bootstrap DGP
    
            for b=1:settings.numboot

                % Generate bootstrap sample based on (possibly lag-augmented) AR estimates
                Y_boot = ar_boot(betahat_ar, res_ar, Y, settings.p+spec_boot_ar.lag_aug, settings.se_homosk, settings.noconst);

                % Estimate on bootstrap sample
                [estims_boot(b,:),ses_boot(b,:)] = ir_estim(Y_boot, rmfield(spec, 'bootstrap'), settings);

            end

        else % Linear regression bootstrap specifications
            
            pseudo_truth = irs;
  
            for h=1:numhorz % Treat each horizon separately

                for b=1:settings.numboot

                    % Generate bootstrap sample
                    [Y_boot, X_boot] = linreg_boot(betahat{h},res{h},X{h},spec.bootstrap,settings.se_homosk);

                    % Run OLS on bootstrap sample
                    [the_linreg_betahat, the_linreg_ses] = linreg(Y_boot,X_boot,settings.se_homosk,true); % Don't add extra intercept
                    estims_boot(b,h) = the_linreg_betahat(1);
                    ses_boot(b,h) = the_linreg_ses(1);

                end

            end
        
        end
        
        % Compute bootstrap confidence intervals
        cis_boot = boot_ci(pseudo_truth', irs', ses', estims_boot, ses_boot, settings.alpha);
        
    end

end