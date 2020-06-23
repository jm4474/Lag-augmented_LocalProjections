function [betahat, se, varcov, res, X_expand] = linreg(Y, X, se_setting, noconst)

    % Linear regression with Eicker-Huber-White s.e.
    
    % Intputs:
    % Y         n x 1   dependent variable data vector
    % X         n x k   covariate data matrix
    % se_setting        EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                   OR function handle: function that returns HAC/HAR sandwich matrix
    % noconst   bool    true: omit intercept
    
    % Outputs:
    % betahat   (k+1) x 1       estimated coefficients
    % se        (k+1) x 1       s.e. of coefficients
    % varcov    (k+1) x (k+1)   var-cov of betahat
    % res       n x 1           residual vector
    % X_expand  n x (k+1)       expanded covariate data matrix with intercept
    
    
    n = size(X,1);
    
    % Determine whether to include intercept
    if noconst
        X_expand = X;
    else
        X_expand = [X ones(n,1)];
    end
    k = size(X_expand,2);
    
    % OLS
    betahat = X_expand\Y;
    
    % Standard errors
    
    if nargout > 1
        
        res = Y-X_expand*betahat;
        XpX = X_expand'*X_expand;
        scores = X_expand.*res;

        if islogical(se_setting)
            if se_setting % If homoskedastic s.e.
                varcov = inv(XpX)*(res'*res)/n;
            else % If EHW s.e.
                varcov = XpX\((scores'*scores)/XpX); % EHW var-cov matrix
            end
        else % If HAC/HAR s.e.
            varcov = XpX\(se_setting(scores)/XpX); % HAC/HAR var-cov matrix
        end

        varcov = n/(n-k)*varcov; % Finite sample adjustment as in Stata
        se = sqrt(diag(varcov)); % Standard errors
    
    end

end