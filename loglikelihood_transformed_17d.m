function loglk=loglikelihood_transformed_17d(X,data)  % Target density being a outer product of 7 univariate densities % Ref: Posterior-guided importance sampling for calculating marginal likelihoods
[r,dim]=size(X);
Y=zeros(r,dim);

Y(:,1)=X(:,1);      % C
Y(:,2)=exp(X(:,2)); % K
Y(:,3)=exp(X(:,3)); % P
Y(:,4)=sqrt(X(:,4).^2+X(:,5).^2); % e
Y(:,5)=mod(atan(X(:,5)./X(:,4))+pi*(X(:,4)<0),2*pi); % omega
Y(:,6)=mod(X(:,6)- Y(:,5),2*pi); % M0
Y(:,7)=exp(X(:,7)); % K
Y(:,8)=exp(X(:,8)); % P
Y(:,9)=sqrt(X(:,9).^2+X(:,10).^2); % e
Y(:,10)=mod(atan(X(:,10)./X(:,9))+pi*(X(:,9)<0),2*pi); % omega
Y(:,11)=mod(X(:,11)- Y(:,10),2*pi); % M0
Y(:,12)=exp(X(:,12)); % K
Y(:,13)=exp(X(:,13)); % P
Y(:,14)=sqrt(X(:,14).^2+X(:,15).^2); % e
Y(:,15)=mod(atan(X(:,15)./X(:,14))+pi*(X(:,14)<0),2*pi); % omega
Y(:,16)=mod(X(:,16)- Y(:,15),2*pi); % M0
Y(:,17)=exp(X(:,17)); % s


loglk= loglikelihood_v5(Y, data);