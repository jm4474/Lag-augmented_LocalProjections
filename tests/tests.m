clear;
addpath('../functions');

% Unit tests

% Run these with Matlab command "runtests" in current directory

% MPM 2020-07-10


%% Impulse responses: AR(2)

lambdas = [0.7 -0.2]; % Eigenvalues of AR polynomial
horzs = 1:10;

A = [sum(lambdas) -prod(lambdas)]; % AR coefficients
irs = var_ir(A,horzs);
c = (lambdas(2)-A(1))/diff(lambdas);
irs2 = c*lambdas(1).^horzs+(1-c)*lambdas(2).^horzs;
assert(norm(irs(:) - irs2(:))<eps);
clear;


%% Impulse responses: VAR(1)

A = [0.9 0; 0.5 0.5];
horzs = 1:10;

irs = var_ir(A,horzs);
n = size(A,1);
nh = length(horzs);
irs2 = zeros(n,n,nh);
for h=1:nh
    irs2(:,:,h) = A^horzs(h);
end
assert(norm(irs(:)-irs2(:))<eps);
clear;


%% Jacobian of impulse responses: VAR(2)

% Compare analytical and numerical derivatives of impulse responses

A = [[0.7 0; 0.5 0.5] [0 0; 0.25 0.25]]; % VAR coefficients
horz = 4; % Horizon of interest
incr = 1e-8; % Increment for numerical derivative

n = size(A,1);
[~, jacob] = var_ir(A,horz);
numparam = numel(A);
the_eye = eye(numparam);
deriv_analy = zeros(n,n,numparam);
deriv_numer = zeros(n,n,numparam);
for i=1:numparam % Loop over parameters in A
    deriv_analy(:,:,i) = reshape(jacob(:,i,1),n,n);
    the_perturb = reshape(the_eye(:,i),size(A));
    deriv_numer(:,:,i) = (var_ir(A+incr*the_perturb,horz)-var_ir(A-incr*the_perturb,horz))/(2*incr); % Two-sided numerical derivative
end
assert(norm(deriv_analy(:)-deriv_numer(:))<1e-7);
clear;


%% Fast Kronecker product

rng(20200710, 'twister');
A = randn(3);
n = 2;
B = randn(size(A,1)*n,2);
assert(norm(kron_fast(A,B,0)-kron(A,eye(n))*B)<1e-10);
assert(norm(kron_fast(A,B,1)-kron(eye(n),A)*B)<1e-10);
clear;


%% Linear regression

rng(20200710, 'twister');
k = 3;
X = randn(100,k);
beta = randn(k,1);
Y = X*beta;
betahat = linreg(Y,X,true,false);
assert(norm(betahat(:)-[beta; 0])<1e-10);
clear;


%% VAR simulation + EWC

T = 1e6;
A = [[0.7 0; 0.5 0.5] [0 0; 0.25 0.25]];
c = [1; 0];
Sigma = [1 0.3; 0.3 1];
rng(20200710, 'twister');

n = size(A,1);
p = size(A,2)/n;
Y = var_sim(A, c, randn(T,n)*chol(Sigma), zeros(p,n));

% Check mean
A1 = eye(n)-sum(reshape(A,n,n,p),3);
assert(norm(mean(Y)'-A1\c)<1e-2);

% Check variance
A_comp = [A; eye(n*(p-1)) zeros(n*(p-1),n)];
Sigma_comp = blkdiag(Sigma, zeros(n*(p-1),n));
var_Y = reshape((eye(n^2*p^2)-kron(A_comp,A_comp))\Sigma_comp(:),n*p,n*p);
assert(norm(cov(Y)-var_Y(1:n,1:n))/norm(var_Y(1:n,1:n))<1e-2);

% Check long-run variance
Omegahat = ewc(Y, round(T^(2/3)))/T;
lrv_Y = (A1\Sigma)/(A1');
assert(norm(Omegahat-lrv_Y)/norm(lrv_Y)<2e-2);
clear;


%% Bias correction: analytical formula

rho = 0.7;
n = 5;
T = 50;

A = rho*eye(n);
Acorr = var_biascorr(A,eye(n),T);
Acorr2 = A + (1+(2+n)*rho)*eye(n)/T;
assert(norm(Acorr-Acorr2)<eps);
clear;


%% Bias correction: simulation

A = [[0.7 0; 0.5 0.5] [0 0; 0.25 0.25]];
Sigma = [1 0.3; 0.3 1];
T = 1e3; % Sample size
numrep = 1e5; % Number of Monte Carlo repetitions
rng(2020710, 'twister');

[n,np] = size(A);
p = np/n;
A_corr = var_biascorr(A, Sigma, T);
b = (A-A_corr)*T; % Analytical negative bias x T

% Simulate data
A_estims = zeros(n,np,numrep);
rngs = randi(2^32-1,numrep,1);
parfor i=1:numrep
    rng(rngs(i), 'twister');
    Y = var_sim(A, [0 0]', mvnrnd([0 0],Sigma,T), zeros(p,n));
    [~,~,betahat_estim] = lp(Y,p-1,1,1:n,true,false);
    A_estims(:,:,i) = betahat_estim(:,1:np);
    if mod(i,ceil(numrep/10))==0
        fprintf('%3d%s\n', 100*i/numrep, '%');
    end
end

assert(norm(b-T*mean(A_estims-A,3))/norm(b)<0.1);
clear;
