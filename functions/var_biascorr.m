function betahat_corr = var_biascorr(betahat, sigma2hat, T)

    % Analytical bias correction for AR(p) estimator
    % Pope (JTSA 1990), equation 9
    
    % Inputs:
    % betahat       p x 1   original AR(p) coefficient estimates
    % sigma2hat     1 x 1   estimate of AR(p) innovation variance
    % T             1 x 1   sample size
    
    % Outputs:
    % betahat_corr  p x 1   bias-corrected coefficient estimates
    
    
    % Set up companion form VAR(1): X_t = A*X_{t-1} + Z_t, where dim(X_t)=p
    p = length(betahat);
    A = [betahat'; eye(p-1), zeros(p-1,1)];
    
    if max(abs(eig(A)))>1 % If original point estimate is outside stationary region, do not bias correct
        betahat_corr = betahat;
        return;
    end
    
    G = diag([sigma2hat, zeros(1,p-1)]); % Var(Z_t)
    Gamma0 = reshape((eye(p^2)-kron(A,A))\G(:),p,p); % Var(X_t)
    
    % Bias correction formula
    aux = inv(eye(p)-A')+(A')/(eye(p)-A'*A');
    lambdas = eig(A);
    for i=1:p
        aux = aux + lambdas(i)*inv(eye(p)-lambdas(i)*A');
    end
    b = G*(aux/Gamma0); % Scaled negative bias
    A_corr = A + b/T; % Bias-corrected VAR(1) coefficients
    
    % If corrected estimate is outside stationary region, reduce bias correction little by little
    % (as recommended by Kilian & Lütkepohl, 2017, ch. 12)
    delta = 1;
    while max(abs(eig(A_corr)))>1 && delta>0
        delta = delta - 0.01;
        A_corr = A + delta*b/T;
    end
    
    % Return bias-corrected AR(p) coefficients
    betahat_corr = A_corr(1,:)';

end