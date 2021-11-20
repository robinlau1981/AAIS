
% multiple thetas, each theta is an numerical or a cell array

function l = loglikelihood_v5(thetas, data)
[r,c]=size(thetas);

% Y(:,1)=thetas(:,1);      % C
% Y(:,2)=exp(thetas(:,2)); % K
% Y(:,3)=exp(thetas(:,3)); % P
% Y(:,4)=sqrt(thetas(:,4).^2+thetas(:,5).^2); % e
% Y(:,5)=mod(atan(thetas(:,5)./thetas(:,4))+pi*(thetas(:,4)<0),2*pi); % omega
% Y(:,6)=mod(thetas(:,6)- Y(:,5),2*pi); % M0
% Y(:,7)=exp(thetas(:,7)); % K
% Y(:,8)=exp(thetas(:,8)); % P
% Y(:,9)=sqrt(thetas(:,9).^2+thetas(:,10).^2); % e
% Y(:,10)=mod(atan(thetas(:,10)./thetas(:,9))+pi*(thetas(:,9)<0),2*pi); % omega
% Y(:,11)=mod(thetas(:,11)- Y(:,10),2*pi); % M0
% Y(:,12)=exp(thetas(:,12)); % K
% Y(:,13)=exp(thetas(:,13)); % P
% Y(:,14)=sqrt(thetas(:,14).^2+thetas(:,15).^2); % e
% Y(:,15)=mod(atan(thetas(:,15)./thetas(:,14))+pi*(thetas(:,14)<0),2*pi); % omega
% Y(:,16)=mod(thetas(:,16)- Y(:,15),2*pi); % M0
% Y(:,17)=exp(thetas(:,17)); % s
% 
% thetas=Y;
l=zeros(r,1);

% for i = 1:r
% 
%   theta = thetas(i,:);
%   if ((((((((theta(2)<0|theta(3)<=0)|theta(6)<0)|theta(6)>=2*pi)|theta(12)<0)|theta(4)<0)|theta(4)>=1)|theta(5)<0)|theta(5)>=2*pi)||((((((((theta(7)<0|theta(8)<=0)|theta(11)<0)|theta(11)>=2*pi)|theta(12)<0)|theta(9)<0)|theta(9)>=1)|theta(10)<0)|theta(10)>=2*pi)
%   
%   % This is meant to yield a very low log likelihood when the parameters are out of the acceptable range.
%     l(i) = -10^100;
%   else
%     l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+theta(12)^2)) - (1/2)*sum(((data.V-model_v4(theta, data.t)).^2)./(data.errors.^2+theta(12)^2));
%     
%   end
% 
% end % for loop

for i = 1:r

  if ((((((((thetas(i,2)<0||thetas(i,3)<=0)||thetas(i,6)<0)||thetas(i,6)>=2*pi)||thetas(i,17)<0)||thetas(i,4)<0)||thetas(i,4)>=1)||thetas(i,5)<0)||thetas(i,5)>=2*pi)||((((((((thetas(i,7)<0||thetas(i,8)<=0)||thetas(i,11)<0)||thetas(i,11)>=2*pi)||thetas(i,17)<0)||thetas(i,9)<0)||thetas(i,9)>=1)||thetas(i,10)<0)||thetas(i,10)>=2*pi)||((((((((thetas(i,12)<0||thetas(i,13)<=0)||thetas(i,16)<0)||thetas(i,16)>=2*pi)||thetas(i,17)<0)||thetas(i,14)<0)||thetas(i,14)>=1)||thetas(i,15)<0)||thetas(i,15)>=2*pi)
  
  % This is meant to yield a very low log likelihood when the parameters are out of the acceptable range.
    l(i) = -Inf;
  else
    l(i) = (-length(data.t)/2)*log(2*pi) - (1/2)*sum(log(data.errors.^2+thetas(i,17)^2)) - (1/2)*sum(((data.V-model_v5(thetas(i,:), data.t)).^2)./(data.errors.^2+thetas(i,17)^2));
    
  end

end % for loop





