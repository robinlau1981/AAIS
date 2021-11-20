
% multiple thetas, each theta is an numerical or a cell array

function l3 = loglikelihood_6d(thetas, data)
[r3,c3]=size(thetas);

% Y(1,j)=X(1,j);      % C
% Y(2,j)=exp(X(2,j)); % K
% Y(3,j)=exp(X(3,j)); % P
% Y(4,j)=sqrt(X(4,j)^2+X(5,j)^2); % e
% Y(5,j)=mod(atan(X(5,j)/X(4,j))+pi*(X(4,j)<0),2*pi); % omega
% Y(6,j)=mod(X(6,j)- Y(5,j),2*pi); % M0
% Y(7,j)=exp(X(7,j)); % s


l3=zeros(r3,1);
for i3 = 1:r3

  if (((((((thetas(i3,2)<0||thetas(i3,3)<=0)||thetas(i3,6)<0)||thetas(i3,6)>=2*pi)||thetas(i3,4)<0)||...
          thetas(i3,4)>=1)||thetas(i3,5)<0)||thetas(i3,5)>=2*pi)
  % This is meant to yield a very low log likelihood when the parameters are out of the acceptable range.
    l3(i3) = -10^100;
  else
    l3(i3) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+10^2)) -...
        (1/2)*sum(((data.V-model_v2(thetas(i3,:), data.t)).^2)./(data.errors.^2+10^2));
  end

end % for loop





