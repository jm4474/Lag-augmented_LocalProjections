clear;
addpath('../functions');

% Test VAR IRF Jacobian formula
% MPM 2020-06-26


%% Settings

A = [[0.7 0; 0.5 0.5] [0 0; 0.25 0.25]]; % VAR coefficients
horz = 4; % Horizon of interest
incr = 1e-8; % Increment for numerical derivative


%% Impulse responses and analytical Jacobian

n = size(A,1);
[~, jacob] = var_ir(A,horz);

numparam = numel(A);
the_eye = eye(numparam);

for i=1:numparam % Loop over parameters in A
    
    fprintf('\n%s%d\n', 'Parameter #', i);
    
    disp('Analytical derivative');
    disp(reshape(jacob(:,i,1),n,n));
    
    disp('Numerical derivative');
    the_perturb = reshape(the_eye(:,i),size(A));
    disp((var_ir(A+incr*the_perturb,horz)-var_ir(A-incr*the_perturb,horz))/(2*incr)); % Two-sided numerical derivative
    
end