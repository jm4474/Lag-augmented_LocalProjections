function [betahat, varcov, res, X_expand] = linreg(Y, X, se_setting, no_const)

    % System linear regression
    % Y_t = beta*X_t + epsilon_t
    % dim(Y_t) = n, dim(beta) = n x k
    
    % Intputs:
    % Y         T x n   dependent variable data matrix
    % X         T x k   covariate data matrix
    % se_setting        EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                   OR function handle: function that returns HAC/HAR sandwich matrix
    % no_const  bool    true: omit intercept
    
    % Outputs:
    % betahat   n x (k+~no_const)       estimated coefficients
    % varcov    (n(k+~no_const)) x      var-cov of vec(betahat)
    %           (n(k+~no_const))                         
    % res       T x n                   residual matrix
    % X_expand  T x (k+~no_const)       covariate data matrix (if no_const)expanded covariate data matrix with intercept
    
    
    [T,n] = size(Y);
    
    % Include intercept if desired
    X_expand = [X ones(T,~no_const)];
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
                varcov = kron_fast(inv(XpX),kron_fast(inv(XpX),scores'*scores,0)',0)'; % EHW var-cov matrix
            end
        else % If HAC/HAR s.e.
            varcov = kron_fast(inv(XpX),kron_fast(inv(XpX),se_setting(scores),0)',0)'; % HAC/HAR var-cov matrix
        end

        varcov = T/(T-k)*varcov; % Finite sample adjustment as in Stata
    
    end

end