function [Y_boot,X_boot] = linreg_boot(beta, res, X, method, homosk)

    % Bootstrap for linear regression
    % Residual-based (i.i.d./wild) or pair/nonparametric bootstrap
    
    % Inputs:
    % beta      k x 1   coefficient estimates in real data
    % res       n x 1   residuals in real data
    % X         n x k   covariate data matrix
    % method    str     'resid': residual bootstrap, 'pair': pair/nonparametric bootstrap
    % homosk    bool    true: homoskedastic bootstrap, false: wild bootstrap (only applies to residual bootstrap)

    
    n = length(res);

    if strcmp(method, 'resid') % Residual bootstrap
        
        if homosk % I.i.d.
            res_boot = res(randi(n,n,1));
        else % Wild
            res_boot = randn(n,1).*res;
        end
        
        Y_boot = X*beta + res_boot;
        X_boot = X;
        
    else % Pair/nonparametric bootstrap
        
        Y = X*beta+res;
        inds = randi(n,n,1);
        Y_boot = Y(inds);
        X_boot = X(inds,:);
        
    end

end