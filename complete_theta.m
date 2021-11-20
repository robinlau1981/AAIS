
% this function takes the given pieces of theta--C, K, P, e, M0, omega-
% and returns a full theta, adding on x=ecos(omga), y=esin(omega), z=M0-omega, and l= log likelihood

function full_theta = complete_theta(theta, data)

full_theta = theta;
full_theta.x = theta.e*cos(theta.omega);
full_theta.y = theta.e*sin(theta.omega);
full_theta.z = mod(theta.omega+theta.M0,2*pi);
full_theta.l = loglikelihood(theta, data);
full_theta.logP = log(theta.P);
full_theta.logK = log(theta.K);
full_theta.logs = log(theta.s);




