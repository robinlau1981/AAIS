function pdf =log_Target_pdf_6d(X,data)  % Target density being a outer product of 7 univariate densities % Ref: Posterior-guided importance sampling for calculating marginal likelihoods
[r,dim]=size(X);
Y=zeros(r,dim);



Y(:,1)=X(:,1);      % C
Y(:,2)=exp(X(:,2)); % K
Y(:,3)=exp(X(:,3)); % P
Y(:,4)=sqrt(X(:,4).^2+X(:,5).^2); % e
Y(:,5)=mod(atan(X(:,5)./X(:,4))+pi*(X(:,4)<0),2*pi); % omega
Y(:,6)=mod(X(:,6)- Y(:,5),2*pi); % M0

loglk= loglikelihood_6d(Y, data);
logprior= log_prior_transformed_6d(X);

pdf =loglk + logprior;
