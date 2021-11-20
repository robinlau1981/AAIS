
% multiple thetas, each theta is an numerical or a cell array

function l = loglikelihood_v3(thetas, data)
[r,c]=size(thetas);

% Y(1,j)=X(1,j);      % C
% Y(2,j)=exp(X(2,j)); % K
% Y(3,j)=exp(X(3,j)); % P
% Y(4,j)=sqrt(X(4,j)^2+X(5,j)^2); % e
% Y(5,j)=mod(atan(X(5,j)/X(4,j))+pi*(X(4,j)<0),2*pi); % omega
% Y(6,j)=mod(X(6,j)- Y(5,j),2*pi); % M0
% Y(7,j)=exp(X(7,j)); % s


l=zeros(r,1);
parfor i = 1:r

  theta = thetas(i,:);

  if ((((((((theta(2)<0||theta(3)<=0)||theta(6)<0)||theta(6)>=2*pi)||theta(7)<0)||theta(4)<0)||theta(4)>=1)||theta(5)<0)||theta(5)>=2*pi)
  % This is meant to yield a very low log likelihood when the parameters are out of the acceptable range.
    l(i) = -Inf;
  else
    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+theta(7)^2)) - (1/2)*sum(((data.V-model_v2(theta, data.t)).^2)./(data.errors.^2+theta(7)^2));
  end

end % for loop





