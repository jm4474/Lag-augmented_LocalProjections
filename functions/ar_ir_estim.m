function [irs, ses, ir_varcov, betahat_estim, res_estim] = ar_ir_estim(Y, p, p_estim, horzs, biascorr, homosk, no_const)
    
    % AR(p) least-squares estimates and delta method s.e.
    % allowing for lag augmentation
    
    % Inputs:
    % Y         T x 1   data vector
    % p         1 x 1   lag length used for impulse response computations
    % p_estim   1 x 1   lag length used for estimation (p_estim >= p)
    % horzs     H x 1   horizons of interest
    % biascorr  bool    true: apply analytical bias correction (Pope, 1990)
    % homosk    bool    true: homoskedastic s.e., false: EHW s.e.
    % no_const  bool    true: omit intercept
    
    % Outputs:
    % irs           H x 1   estimated impulse responses at select horizons
    % ses           H x 1   s.e. of impulse responses at select horizons
    % ir_varcov     H x H   var-cov matrix of impulse responses at select horizons
    % betahat_estim p_estim x 1     coefficient estimates (possibly bias-corrected, possibly including extra intercept)
    % res_estim     T-p_estim x 1   estimation residuals
    
    
    % One-step forecasting regression of Y_{t+1} on (Y_t, ..., Y_{t-p_estim+1})
    [~,~,betahat_estim,betahat_estim_varcov,res_estim] = lp(Y,p_estim-1,1,homosk,no_const);
    
    % If bias correction is desired...
    if biascorr
        sigma2hat = (res_estim'*res_estim)/(length(res_estim)-p_estim-1); % Residual variance estimate
        betahat_estim(1:p_estim) = ar_biascorr(betahat_estim(1:p_estim), sigma2hat, length(Y));
    end
    
    % Only use first p AR coefficients to compute impulse responses
    betahat = betahat_estim(1:p);
    betahat_varcov = betahat_estim_varcov(1:p,1:p);
    
    if nargout==1
        irs = ar_ir(betahat,horzs); % Compute impulse responses
    else
        [irs, jacob] = ar_ir(betahat,horzs); % Compute impulse responses and Jacobian
        ir_varcov = jacob*betahat_varcov*jacob';
        ses = sqrt(diag(ir_varcov)); % S.e.
    end

end