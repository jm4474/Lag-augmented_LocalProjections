function cis = boot_ci(pseudo_truth, estim_data, ses_data, estims_boot, ses_boot, alpha)

    % Bootstrap confidence intervals
    
    % Inputs:
    % pseudo_truth  1 x k   pseudo-true parameter vector used to generate bootstrap samples
    % estim_data    1 x k   parameter estimates in real data
    % ses_data      1 x k   s.e. in real data
    % estims_boot   B x k   bootstrapped parameter estimates
    % ses_boot      B x k   bootstrapped s.e.
    % alpha         1 x 1   significance level
    
    % Outputs:
    % cis           2 x k x 3   bootstrap intervals (1st index: lower and upper limit, 2nd index: parameter, 3rd index: type of interval - Efron, Hall, or Hall percentile-t)
    
    
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
    cis = reshape([ci_efron ci_hall ci_hall_t], 2, [], 3);

end