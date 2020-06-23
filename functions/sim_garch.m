function U = sim_garch(alpha_0,alpha_1,beta_1,T)

    % Simulate GARCH(1,1) time series
    
    % Model:
    % u_t = sigma_t*epsilon_t
    % sigma_t^2 = alpha_0 + alpha_1*u_{t-1}^2 + beta_1*sigma_{t-1}^2
    %           = alpha_0 + (alpha_1*epsilon_{t-1}^2 + beta_1)*sigma_{t-1}^2
    
    % Inputs:
    % alpha_0   1 x 1   intercept parameter
    % alpha_1   1 x 1   MA parameter
    % beta_1    1 x 1   AR parameter
    % T         1 x 1   sample size
    
    % Outputs
    % U         T x 1   simulated GARCH process u_t
    
    
    epsilons = randn(T,1); % Shocks
    
    sigma2s = zeros(T,1);
    sigma2s(1) = alpha_0/(1-alpha_1-beta_1); % Initialize sigma_t^2 at its ergodic mean
    
    for t=2:T
        sigma2s(t) = alpha_0 + (alpha_1*epsilons(t-1)^2 + beta_1)*sigma2s(t-1);
    end
    
    U = sqrt(sigma2s).*epsilons;

end