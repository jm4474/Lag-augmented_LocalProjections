function [Y_boot,X_boot] = linreg_boot(beta, res, X, method, homosk)

    % Bootstrap for linear regression
    % Residual-based (i.i.d./wild) or pair/nonparametric bootstrap
    
    % Inputs:
    % beta      k x 1   coefficient estimates in real data
    % res       T x 1   residuals in real data
    % X         T x k   covariate data matrix
    % method    str     'resid': residual bootstrap, 'pair': pair/nonparametric bootstrap
    % homosk    bool    true: homoskedastic bootstrap, false: wild bootstrap (only applies to residual bootstrap)

    
    T = length(res);

    if strcmp(method, 'resid') % Residual bootstrap
        
        if homosk % I.i.d.
            res_boot = res(randi(T,T,1));
        else % Wild
            res_boot = randn(T,1).*res;
        end
        
        Y_boot = X*beta + res_boot;
        
        X_boot = X;
        
    else % Pair/nonparametric bootstrap
        
        Y      = X*beta + res;
        
        inds   = randi(T,T,1);
        
        Y_boot = Y(inds);
        
        X_boot = X(inds,:);
        
    end

end