clear;
addpath('../functions');

% Test bias correction function
% MPM 2020-06-26


%% Settings

% DGP
A = [[0.7 0; 0.5 0.5] [0 0; 0.25 0.25]];
Sigma = [1 0.3; 0.3 1];

T = 1e3; % Sample size
numrep = 1e5; % Number of Monte Carlo repetitions
rng(202006261);


%% Analytical bias

[n,np] = size(A);
p = np/n;

A_corr = var_biascorr(A, Sigma, T);
b = (A-A_corr)*T; % Analytical bias x T


%% Simulate bias

A_estims = zeros(n,np,numrep);

parfor i=1:numrep
    
    Y = var_sim(A, [0 0]', mvnrnd([0 0],Sigma,T), zeros(p,n));
    [~,~,betahat_estim] = lp(Y,p-1,1,1:n,true,false);
    A_estims(:,:,i) = betahat_estim(:,1:np);
    
    if mod(i,ceil(numrep/10))==0
        fprintf('%3d%s\n', 100*i/numrep, '%');
    end
    
end


%% Compare analytical and simulation bias

disp('Analytical bias');
disp(b);
disp('Simulated bias');
disp(T*mean(A_estims-A,3));

