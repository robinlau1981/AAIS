% computes the log likelihood of data, parameters

% I'm now modifying this function so that it can take a theta_list--a list containing theta structs--and return a list of loglikelihoods for all of them.

function l = loglikelihood_0p_v2(thetas, data)

[r,c]=size(thetas);
l=zeros(r,1);
for i = 1:r

    

    %  if (theta(i,2)<0)
    %    l(i) = -10^100;
    %  else
    %    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+theta(i,2)^2)) - (1/2)*sum(((data.V-theta(i,1)).^2)./(data.errors.^2+theta(i,2)^2));
    %  end

    %theta(i,2) = exp(theta.logs);

    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+thetas(i,2)^2)) - (1/2)*sum(((data.V-thetas(i,1)).^2)./(data.errors.^2+thetas(i,2)^2));

end % for loop