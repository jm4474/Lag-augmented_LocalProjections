function [ir, ir_varcov, betahat, betahat_varcov, res, X] = lp(Y,num_lags,horz,resp_ind,se_setting,no_const)

    % Local projection with lagged controls
    
    % Inputs:
    % Y         T x n       data matrix
    % num_lags  1 x 1       number of lags of Y to control for (in addition
    % to the contemporaneous regressor); in comments refered to as p
    % horz      1 x 1       horizon of interest
    % resp_ind  1 x m       indices of response variables Y_t of interest
    % se_setting            EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                       OR function handle: function that returns HAC/HAR sandwich matrix
    % no_const   bool       true: omit intercept
    
    % Outputs:
    % ir                m x n               estimated impulse responses at select horizons
    % ir_varcov         mn x mn             var-cov of vec(ir)
    % betahat           (np+n+1) x 1        full vector of estimated regression coefficients (+1 if
    % no_const is false)
    % betahat_varcov    (np+n+1) x (np+n+1) var-cov of betahat (+1 if
    % no_const is false)
    % res               (T-p) x m           residuals
    % X                 (T-p) x (np+n+1)    covariate data matrix including intercept
    
    
    % Covariate matrix
    Y_lag = lagmatrix(Y,0:num_lags);
    X = Y_lag(num_lags+1:end-horz,:);
    
    % Local projection
    [betahat, betahat_varcov, res, X] = linreg(Y(num_lags+horz+1:end,resp_ind),X,se_setting,no_const);
    
    n = size(Y,2);
    m = length(resp_ind);
    ir = betahat(:,1:n);
    ir_varcov = betahat_varcov(1:m*n,1:m*n);

end