function A_corr = var_biascorr(A, Sigma, T)

    % Analytical bias correction for VAR(p) estimator
    % Pope (JTSA 1990), equation 9
    
    % Inputs:
    % A             n x np  original VAR(p) coefficient estimates [A_1,...,A_p]
    % Sigma     	n x n   estimate of VAR(p) innovation variance
    % T             1 x 1   sample size
    
    % Outputs:
    % A_corr        n x np  bias-corrected coefficient estimates
    
    
    % Set up companion form: X_t = A*X_{t-1} + Z_t, where dim(X_t)=n*p
    [n,np] = size(A);
    A_comp = [A; eye(np-n), zeros(np-n,n)];
    
    if max(abs(eig(A_comp)))>1 % If original point estimate is outside stationary region, do not bias correct
        A_corr = A;
        return;
    end
    
    G = blkdiag(Sigma, zeros(np-n)); % Var(Z_t)
    Gamma0 = dlyap(A_comp, G); % Var(X_t), requires Control System Toolbox
    % The following slower command does not require Control System Toolbox:
%     Gamma0 = reshape((eye(np^2)-kron(A_comp,A_comp))\G(:),np,np);  
    
    % Bias correction formula
    aux = inv(eye(np)-A_comp')+(A_comp')/(eye(np)-A_comp'*A_comp');
    lambdas = eig(A_comp);
    for lamb = lambdas'
        aux = aux + lamb*inv(eye(np)-lamb*A_comp');
    end
    b = G*(aux/Gamma0); % Scaled negative bias
    A_corr = A_comp + b/T; % Bias-corrected companion form coefficients
    
    % If corrected estimate is outside stationary region, reduce bias correction little by little
    % (as recommended by Kilian & Lütkepohl, 2017, ch. 12)
    delta = 1;
    while max(abs(eig(A_corr)))>1 && delta>0
        delta = delta - 0.01;
        A_corr = A_comp + delta*b/T;
    end
    
    % Return bias-corrected VAR(p) coefficients
    A_corr = A_corr(1:n,:);

end