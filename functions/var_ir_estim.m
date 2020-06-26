function [irs, irs_varcov, Ahat_estim, res_estim] = var_ir_estim(Y, p, p_estim, horzs, bias_corr, homosk, no_const)
    
    % VAR(p) least-squares estimates and delta method s.e.
    % allowing for lag augmentation
    
    % Inputs:
    % Y         T x n   data vector
    % p         1 x 1   lag length used for impulse response computations
    % p_estim   1 x 1   lag length used for estimation (p_estim >= p)
    % horzs     H x 1   horizons of interest
    % bias_corr bool    true: apply analytical bias correction (Pope, 1990)
    % homosk    bool    true: homoskedastic s.e., false: EHW s.e.
    % no_const  bool    true: omit intercept
    
    % Outputs:
    % irs           n x n x H           estimated impulse responses Theta_h at select horizons
    % irs_varcov    n^2 x n^2 x H       var-cov matrices of vec(Theta_h) at select horizons
    % Ahat_estim    n x np              VAR coefficient estimates [A_1,...,A_p] (possibly bias-corrected, possibly including intercept as last column)
    % res_estim     (T-p_estim) x n     estimation residuals
    
    
    [T,n] = size(Y);
    
    % One-step forecasting regression of Y_{t+1} on (Y_t, ..., Y_{t-p_estim+1})
    [~,~,Ahat_estim,Ahat_estim_varcov,res_estim] = lp(Y,p_estim-1,1,1:n,homosk,no_const);
    
    % If bias correction is desired...
    if bias_corr
        Sigmahat = (res_estim'*res_estim)/(size(res_estim,1)-n*p_estim-1+no_const); % Residual variance estimate
        Ahat_estim(:,1:end-1+no_const) = var_biascorr(Ahat_estim(:,1:end-1+no_const), Sigmahat, T);
    end
    
    % Only use first p VAR coefficient matrices to compute impulse responses
    Ahat = Ahat_estim(:,1:n*p);
    Ahat_varcov = Ahat_estim_varcov(1:n^2*p,1:n^2*p);
    
    if nargout==1
        irs = var_ir(Ahat,horzs); % Compute impulse responses
    else
        [irs, jacob] = var_ir(Ahat,horzs); % Compute impulse responses and Jacobian
        nh = length(horzs);
        irs_varcov = zeros(n^2,n^2,nh);
        for h=1:nh
            irs_varcov(:,:,h) = jacob(:,:,h)*Ahat_varcov*jacob(:,:,h)'; % Var-cov for impulse response matrix at this horizon
        end
    end

end