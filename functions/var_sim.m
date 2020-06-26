function Y = var_sim(A, c, U, Y_init)

    % Simulate VAR(p) data
    % given innovations and initial conditions
    
    % Inputs:
    % A         n x np  VAR coefficients [A_1,...,A_p]
    % c         n x 1   Intercepts
    % U         T x n   Innovations data matrix
    % Y_init    p x n   Initial conditions [y_1,...,y_p]' (only last p rows will be used)
    
    % Outputs:
    % Y         T x n   Simulated data matrix [y_1,...,y_T]', discarding initial obs.
    
    
    % Dimensions
    [T,n] = size(U);
    p = size(A,2)/n;
    
    % Iterate
    Y = zeros(T+p,n);
    Y(1:p,:) = Y_init(end-p+1:end,:); % Initial conditions
    for t=p+1:T+p
        Y(t,:) = c' + reshape(Y(t-1:-1:t-p,:)',1,n*p)*A' + U(t-p,:);
    end
    Y = Y(p+1:end,:); % Discard initial conditions

end