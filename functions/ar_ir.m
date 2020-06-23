function [irs, jacobs] = ar_ir(beta,horzs)

    % AR(p) impulse responses and Jacobian wrt. parameters
    
    % Inputs:
    % beta      p x 1   AR coefficients (beta_1, ..., beta_p)
    % horzs     H x 1   horizons of interest
    
    % Outputs:
    % irs       H x 1   impulse responses at select horizons
    % jacobs    H x p   Jacobian of impulse responses at select horizons wrt. beta components
    
    
    % Dimensions
    nh = length(horzs);
    maxh = max(horzs);
    p = length(beta);
    
    irs = zeros(nh,1); % Will contain impulse responses at select horizons
    jacobs = zeros(nh,p); % Will contain Jacobian of impulse responses at select horizons wrt. AR coefficients
    
    ir_p = [1; zeros(p-1,1)]; % Will contain last p impulse responses
    jacob_p = zeros(p); % Will contain Jacobians of the last p impulse responses wrt. AR coefficients
    
    for i=1:maxh % Loop through horizons
        
        the_beta = beta(1:min(i,p));
        the_past_ir = ir_p(1:min(i,p));
        the_ir = the_beta'*the_past_ir; % Impulse response at horizon i
        ir_p = [the_ir; ir_p(1:end-1)]; % Shift forward in time
        
        irs(horzs==i) = the_ir; % Store impulse response at select horizons
        
        if nargout>1
            the_jacob_p = zeros(1,p);
            the_jacob_p(1:min(i,p)) = the_past_ir' + the_beta'*jacob_p(1:min(i,p),1:min(i,p));
            jacob_p = [the_jacob_p; jacob_p(1:end-1,:)]; % Shift forward in time
            
            the_ind = find(horzs==i,1);
            if ~isempty(the_ind)
                jacobs(the_ind,:) = the_jacob_p; % Store Jacobian at select horizons
            end
        end
        
    end

end