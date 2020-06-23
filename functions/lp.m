function [ir, se, betahat, betahat_varcov, res, X] = lp(Y,num_lags,horz,se_setting,noconst)

    % Local projection with lagged controls
    
    % Inputs:
    % Y         T x 1       data vector
    % num_lags  1 x 1       number of lags of Y to control for (in addition to the contemporaneous regressor)
    % horz      1 x 1       horizon of interest
    % se_setting        EITHER bool: if true, homoskedastic s.e.; if false, EHW s.e.
    %                   OR function handle: function that returns HAC/HAR sandwich matrix
    % noconst   bool    true: omit intercept
    
    % Outputs:
    % ir        1 x 1               estimated impulse response at select horizon
    % se        1 x 1               s.e. of impulse response at select horizon
    % betahat   (p+2) x 1           full vector of estimated regression coefficients
    % betahat_varcov (p+2) x (p+2)  var-cov of betahat
    % res       (T-p) x 1           residuals
    % X         (T-p) x (p+2)       covariate data matrix including intercept
    
    
    % Covariate matrix
    Y_lag = lagmatrix(Y,0:num_lags);
    X = Y_lag(num_lags+1:end-horz,:);
    
    % Local projection
    [betahat, ses, betahat_varcov, res, X] = linreg(Y(num_lags+horz+1:end),X,se_setting,noconst);
    ir = betahat(1);
    se = ses(1);

end