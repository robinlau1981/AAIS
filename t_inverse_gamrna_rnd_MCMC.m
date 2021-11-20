% Generate samples from inverse gamma distribution using MCMC

clear all

% number of iterations
n = 10000;
% burn-in rounds
burnin = 2000;

vari = randn(1, n);

%str = 'beta^alpha/gamma(alpha)*(1/x)^(alpha + 1)*exp(-beta/x)'; % inverse-gamma distribution
%inverGam = inline(str, 'x', 'alpha', 'beta');

str = 'beta^alpha/gamma(1)*(1/x)^(alpha + 1)*exp(-beta/x)'; % inverse-gamma distribution
inverGam = inline(str, 'x', 'alpha', 'beta');
alpha = 1;
beta = 10;

str = 'exp(-0.5*((x-mu)/sig).^2)';
norm = inline(str, 'x', 'mu', 'sig');
%str = 'exp(-0.5*((x-10)/sig).^2)';
sig = 1;

vari(1) = .1; % initial value
for index = 2:n
    % normal distribution is used as proposal distribution
    y = vari(index - 1) + randn(1)*sig;
    % draw u from normal distribution
    u = rand;
    % calculate alpha
    % Metropolis-Hasting algorithm
    if y>0
    r = min(1,inverGam(y',alpha, beta)*normpdf(vari(index - 1), y,sig)/...
        (inverGam(vari(index-1)', alpha, beta))*normpdf(y,vari(index - 1), sig));
    else
        r=0;
    end
    if u <= r
        vari(index) = y;
    else
        vari(index) = vari(index - 1);
    end
end

% remove burn-in elements
vari = vari(burnin:n);
hist(vari,1e2);