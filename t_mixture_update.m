function [W_m_minus,Mu_minus,Sigma_minus]=t_mixture_update(M,W_m,Mu,Sigma,X,q,df,rp,proposal)
[N,dim]=size(X);
W_m_minus=zeros(1,M);
Mu_minus=zeros(M,dim);
Sigma_minus=zeros(dim,dim,M);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
rou=zeros(M,N);
u=zeros(M,N);
parfor ii=1:M
    %rou(ii,:)=W_m(ii)*((t_pdf(X,Mu(ii,:),Sigma(:,:,ii),df))./proposal)';
    rou(ii,:)=W_m(ii)*(rp(:,ii)./proposal)';
    W_m_minus(ii)=rou(ii,:)*q;
    u(ii,:)=(df+dim)*ones(1,N)./(df*ones(1,N)+MahalanobisDist(X,Mu(ii,:),Sigma(:,:,ii)));
    Mu_minus(ii,:)=rou(ii,:).*u(ii,:)*(repmat(q,1,dim).*X)/(rou(ii,:).*u(ii,:)*q);%W_m_minus(ii);

    %Sigma_minus(:,:,ii)=(repmat(q,1,dim).*(X-repmat(Mu_minus(ii,:),N,1)))'*((X-repmat(Mu_minus(ii,:),N,1)).*repmat((rou(ii,:).*u(ii,:))',1,dim))/W_m_minus(ii);
    Sigma_minus(:,:,ii)=(X-repmat(Mu_minus(ii,:),N,1))'*((X-repmat(Mu_minus(ii,:),N,1)).*repmat(q,1,dim).*repmat(rou(ii,:)',1,dim).*repmat(u(ii,:)',1,dim))/W_m_minus(ii);
    %Sigma_minus(:,:,ii)=diag(sum(repmat(q,1,dim).*((X-repmat(Mu_minus(ii,:),N,1)).^2).*repmat((rou(ii,:)'.*u(ii,:)'),1,dim))/W_m_minus(ii));
    %         for i=1:N
    %             Sigma_minus(:,:,ii)=Sigma_minus(:,:,ii)+q(i,j)*(X(i,:)-Mu_minus(ii,:))'*(X(i,:)-Mu_minus(ii,:))*rou(ii,i)*u(ii,i)/W_m_minus(ii);
    %         end
    %         parfor d=1:dim
    %             Sigma_minus(d,d,ii)=(sum(q'.*((X(:,d)-Mu_minus(ii,d)).^2)'.*(rou(ii,:).*u(ii,:)))+2*beta(d)/N)/(W_m_minus(ii)+2*(alpha+1)/N);
    %             %Sigma_minus(d,d,ii)=(q.*(X(:,d)-Mu_minus(ii,d)))'*((X(:,d)-Mu_minus(ii,d)).*(rou(ii,:).*u(ii,:))')/W_m_minus(ii);
    %         end
end