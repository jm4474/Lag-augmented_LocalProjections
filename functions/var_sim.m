function Y = var_sim(beta, c, U, Y_init)

    % Simulate VAR(p) data
    % given innovations and initial conditions
    
    
    % Dimensions
    [T,n] = size(U);
    p = size(beta,2)/n;
    
    % Iterate
    Y = zeros(T+p,n);
    Y(1:p,:) = Y_init(end-p+1:end,:); % Initial conditions
    for t=p+1:T+p
        Y(t,:) = c' + reshape(Y(t-1:-1:t-p,:)',1,n*p)*beta' + U(t-p,:);
    end
    Y = Y(p+1:end,:); % Discard initial conditions

end