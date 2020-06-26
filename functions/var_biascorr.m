function betahat_corr = var_biascorr(betahat, Sigmahat, T)

    % Analytical bias correction for VAR(p) estimator
    % Pope (JTSA 1990), equation 9
    
    % Inputs:
    % betahat       n x np  original VAR(p) coefficient estimates [A_1,...,A_p]
    % Sigmahat     	n x n   estimate of VAR(p) innovation variance
    % T             1 x 1   sample size
    
    % Outputs:
    % betahat_corr  n x np  bias-corrected coefficient estimates
    
    
    % Set up companion form: X_t = A*X_{t-1} + Z_t, where dim(X_t)=n*p
    [n,np] = size(betahat);
    A = [betahat; eye(np-n), zeros(np-n,n)];
    
    if max(abs(eig(A)))>1 % If original point estimate is outside stationary region, do not bias correct
        betahat_corr = betahat;
        return;
    end
    
    G = blkdiag(Sigmahat, zeros(np-n)); % Var(Z_t)
    Gamma0 = reshape((eye(np^2)-kron(A,A))\G(:),np,np); % Var(X_t)
    
    % Bias correction formula
    aux = inv(eye(np)-A')+(A')/(eye(np)-A'*A');
    lambdas = eig(A);
    for lamb = lambdas'
        aux = aux + lamb*inv(eye(np)-lamb*A');
    end
    b = G*(aux/Gamma0); % Scaled negative bias
    A_corr = A + b/T; % Bias-corrected companion form coefficients
    
    % If corrected estimate is outside stationary region, reduce bias correction little by little
    % (as recommended by Kilian & Lütkepohl, 2017, ch. 12)
    delta = 1;
    while max(abs(eig(A_corr)))>1 && delta>0
        delta = delta - 0.01;
        A_corr = A + delta*b/T;
    end
    
    % Return bias-corrected VAR(p) coefficients
    betahat_corr = A_corr(1:n,:);

end