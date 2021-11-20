function [W_m_minus,Mu_minus,Sigma_minus,df]=t_mixture_update_v2(M,W_m,Mu,Sigma,X,q,df,rp,proposal,tao)
[N,dim]=size(X);
W_m_minus=zeros(1,M);
Mu_minus=zeros(M,dim);
Sigma_minus=zeros(dim,dim,M);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%rou=zeros(M,N);
%u=zeros(M,N);
for ii=1:M    
    rou=W_m(ii)*(rp(:,ii)./proposal)';
    W_m_minus(ii)=rou*q;
    u=(df(ii)+dim)*ones(1,N)./(df(ii)*ones(1,N)+MahalanobisDist(X,Mu(ii,:),Sigma(:,:,ii)));
    Mu_minus(ii,:)=rou.*u*(repmat(q,1,dim).*X)/(rou.*u*q);%W_m_minus(ii);

    Sigma_minus(:,:,ii)=(X-repmat(Mu_minus(ii,:),N,1))'*((X-repmat(Mu_minus(ii,:),N,1)).*...
        repmat(q,1,dim).*repmat(rou',1,dim).*repmat(u',1,dim))/(rou.*u*q);%W_m_minus(ii);
    
%      for d=1:dim
%                 %Sigma_minus(d,d,ii)=(sum(q'.*((X(:,d)-Mu_minus(ii,d)).^2)'.*(rou(ii,:).*u(ii,:)))+2*beta(d)/N)/(W_m_minus(ii)+2*(alpha+1)/N);
%                 Sigma_minus(d,d,ii)=(sum(q'.*((X(:,d)-Mu_minus(ii,d)).^2)'.*(rou.*u)))/(W_m_minus(ii));
%                 %Sigma_minus(d,d,ii)=(q.*(X(:,d)-Mu_minus(ii,d)))'*((X(:,d)-Mu_minus(ii,d)).*(rou(ii,:).*u(ii,:))')/W_m_minus(ii);
%      end
    
    %N_tmp=W_m_minus(ii)*N;
    df(ii)=df(ii);
    %df(ii)=fzero(@(v) -psi(.5*v)+log(.5*v)+1+1/N_tmp*(rou.*(log(u)-tao'))*q...
    %   +psi(.5*(v+dim))-log(.5*(v+dim)),df(ii)); 
    % based on 'Robust mixture modelling using the t distribution'
end