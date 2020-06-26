function Y_boot = var_boot(A, res, Y, p, homosk, no_const)

    % VAR residual bootstrap, homoskedastic or wild
    
    % Inputs:
    % A         n x np      VAR(p) coefficient matrices [A_1,...,A_p]
    %                       OPTIONAL: column dimension could exceed n*p, in which case last column equals intercept
    %                       NOTE: only A(:,1:n*p) and A(:,end) will be used
    % res       T_res x n   residuals
    % Y         T x n       data vector
    % p         1 x 1       lag length
    % homosk    bool        true: homoskedastic bootstrap, false: wild bootstrap
    % no_const  bool        true: exclude intercept from bootstrap samples
    
    % Outputs:
    % Y_boot    T x n       bootstrap sample
    

    % Dimensions
    [T,n] = size(Y);
    T_res = size(res,1); % Effective sample size for residuals (need not equal T-p)
    
    % Draw block of initial T-T_res obs. from real data
    ind_init = randi(T_res+1);
    Y_init = Y(ind_init:ind_init+T-T_res-1,:);
    
    % Intercept (if applicable)
    c = zeros(n,1);
    if ~no_const
        c = A(:,end);
    end

    % Draw residuals
    if homosk % I.i.d. bootstrap
        res_boot = res(randi(T_res,T_res,1),:);
    else % Wild bootstrap
        res_boot = randn(T_res,1).*res;
    end
    
    % Generate VAR(p) data, with residuals and initial conditions as above
    Y_boot = [Y_init; var_sim(A(:,1:n*p), c, res_boot, Y_init)];

end