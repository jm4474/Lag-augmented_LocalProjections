function [irs, jacob] = var_ir(A,horzs)

    % VAR(p) reduced-form impulse responses and Jacobian wrt. parameters
    
    % Inputs:
    % A         n x np  VAR coefficient matrices (A_1, ..., A_p)
    % horzs     1 x H   horizons of interest
    
    % Outputs:
    % irs       n x n x H           reduced-form impulse responses Theta_h at select horizons
    % jacob     n^2 x (n^2*p) x H   Jacobian of vec(Theta_h) at select horizons wrt. vec(A)
    
    
    % Dimensions
    nh = length(horzs);
    maxh = max(horzs);
    [n,np] = size(A);
    p = np/n;
    
    irs = zeros(n,n,nh); % Will contain impulse responses at select horizons
    jacob = zeros(n^2,n^2*p,nh); % Will contain Jacobian of impulse responses at select horizons wrt. vec(A)
    
    ir_p = [eye(n); zeros(n*(p-1),n)]; % Will contain last p impulse responses, stacked vertically
    jacob_p = zeros(n^2,n^2*p,p); % Will contain last p values of the Jacobian of vec(Theta_h) wrt. vec(A)
    
    for h=1:maxh % Loop through horizons
        
        the_A = A(:,1:n*min(h,p));
        the_past_ir = ir_p(1:n*min(h,p),:);
        the_ir = the_A*the_past_ir; % Impulse response at horizon h
        ir_p = [the_ir; ir_p(1:end-n,:)]; % Shift forward in time
        
        the_ind = find(horzs==h,1);
        if ~isempty(the_ind)
            irs(:,:,the_ind) = the_ir; % Store impulse response at select horizons
        end
        
        if nargout>1
            
            % Jacobian at horizon h via chain rule
            the_jacob_p = zeros(n^2,n^2*p);
            the_jacob_p(:,1:n^2*min(h,p)) = kron(the_past_ir',eye(n));
            for l=1:min(h,p)
                the_jacob_p(:,1:n^2*min(h,p)) = the_jacob_p(:,1:n^2*min(h,p)) + kron_fast(A(:,(l-1)*n+1:l*n),jacob_p(:,1:n^2*min(h,p),l),1);
            end
            
            % Shift forward in time
            jacob_p(:,:,2:end) = jacob_p(:,:,1:end-1);
            jacob_p(:,:,1) = the_jacob_p;
            
            if ~isempty(the_ind)
                jacob(:,:,the_ind) = the_jacob_p; % Store Jacobian at select horizons
            end
            
        end
        
    end

end