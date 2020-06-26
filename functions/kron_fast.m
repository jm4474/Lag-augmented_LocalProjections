function C = kron_fast(A,B,type)

    % Efficient computation of
    % C = kron(A, eye(n))*B (type=0)
    % or
    % C = kron(eye(n), A)*B (type=1)
    
    % From "MATLAB array manipulation tips and tricks"
    % by Peter J. Acklam (2002), section 10.8
    
    
    % Matrix dimensions
    [p,q] = size(A);
    [qn,m] = size(B);
    n = qn/q;
    
    % Compute
    switch type
        case 0
            C = reshape(reshape(B.', [n*m q])*A.', [m p*n]).';
        case 1
            C = reshape(A*reshape(B, [q n*m]), [p*n m]);
    end

end