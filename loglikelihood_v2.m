
% computes the log likelihood of data, parameters

% I'm now modifying this function so that it can take a theta_list--a list containing theta structs--and return a list of loglikelihoods for all of them.

% data is a struct with elements t, V, and e (time, radial velocity, and
% errors.  All are vectors.)

function l = loglikelihood_v2(thetas, data)

for i = 1:length(thetas)

  theta = thetas(i);

  if ((((((((theta.K<0|theta.P<=0)|theta.M0<0)|theta.M0>=2*pi)|theta.s<0)|theta.e<0)|theta.e>=1)|theta.omega<0)|theta.omega>=2*pi)|((((((((theta.K2<0|theta.P2<=0)|theta.M02<0)|theta.M02>=2*pi)|theta.s<0)|theta.e2<0)|theta.e2>=1)|theta.omega2<0)|theta.omega2>=2*pi)
  % This is meant to yield a very low log likelihood when the parameters are out of the acceptable range.
    l(i) = -10^100;
  else
    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+theta.s^2)) - (1/2)*sum(((data.V-model_v3(theta, data.t)).^2)./(data.errors.^2+theta.s^2));
  end

end % for loop





