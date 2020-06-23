function Y_boot = ar_boot(beta, res, Y, p, homosk, noconst)

    % AR residual bootstrap, homoskedastic or wild
    
    % Inputs:
    % beta      m x 1       AR(m) coefficient estimates
    %                       OPTIONAL: length could exceed p, in which case last element equals intercept
    %                       NOTE: only beta(1:p) and beta(end) will be used
    % res       T_res x 1   residuals
    % Y         T x 1       data vector
    % p         1 x 1       lag length
    % homosk    bool        true: homoskedastic bootstrap, false: wild bootstrap
    % noconst   bool        true: exclude intercept from bootstrap samples
    
    % Outputs:
    % Y_boot    T x 1       bootstrap sample
    

    % Dimensions
    T = length(Y); % Sample size
    T_res = length(res); % Effective sample size for residuals (need not equal T-p)
    
    % Draw block of initial T-T_res obs. from real data
    ind_init = randi(T_res+1);
    Y_init = Y(ind_init:ind_init+T-T_res-1);
    
    % Intercept (if applicable)
    nu = 0;
    if ~noconst
        nu = beta(end);
    end

    % Draw residuals
    if homosk % I.i.d. bootstrap
        res_boot = res(randi(T_res,T_res,1));
    else % Wild bootstrap
        res_boot = randn(T_res,1).*res;
    end
    
    % Generate AR(p) data, ensuring correct initial conditions (drawn above)
    aux = triu(toeplitz(Y_init(end-p+1:end)));
    Y_boot = filter(1, [1 -beta(1:p)'], [Y_init(end-p+1:end) - [zeros(p,1), aux(:,1:end-1)]'*beta(1:p); nu + res_boot]);
    Y_boot = [Y_init(1:T-T_res-p); Y_boot];

end