% computes the log likelihood of data, parameters

% I'm now modifying this function so that it can take a theta_list--a list containing theta structs--and return a list of loglikelihoods for all of them.

function l = loglikelihood_0p(thetas, data)

for i = 1:length(thetas)

    theta = thetas(i);

    %  if (theta.s<0)
    %    l(i) = -10^100;
    %  else
    %    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+theta.s^2)) - (1/2)*sum(((data.V-theta.C).^2)./(data.errors.^2+theta.s^2));
    %  end

    %theta.s = exp(theta.logs);

    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+theta.s^2)) - (1/2)*sum(((data.V-theta.C).^2)./(data.errors.^2+theta.s^2));

end % for loop