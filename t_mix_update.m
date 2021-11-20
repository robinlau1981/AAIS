function MixtureUpdated=t_mix_update(X,Proposal,df)

M=Proposal.M;
W_m=Proposal.W;
Mu=Proposal.Mu;
Sigma=Proposal.Sigma;
q=X.NormalizedWeight;
rp=X.Resp;
proposal=X.Proposal;

[N,dim]=size(X.Values);
W_m_minus=zeros(1,M);
Mu_minus=zeros(M,dim);
Sigma_minus=zeros(dim,dim,M);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
rou=zeros(M,N);
u=zeros(M,N);
for ii=1:M    
    rou(ii,:)=W_m(ii)*(rp(:,ii)./proposal)';
    W_m_minus(ii)=rou(ii,:)*q;    
    u(ii,:)=(df+dim)*ones(1,N)./(df*ones(1,N)+sqdist(X.Values',Mu(ii,:)',inv_posdef(Sigma(:,:,ii)))');
    Mu_minus(ii,:)=rou(ii,:).*u(ii,:)*(repmat(q,1,dim).*X.Values)/(rou(ii,:).*u(ii,:)*q);%W_m_minus(ii);
    Sigma_minus(:,:,ii)=(X.Values-repmat(Mu_minus(ii,:),N,1))'*((X.Values-repmat(Mu_minus(ii,:),N,1)).*repmat(q,1,dim).*repmat(rou(ii,:)',1,dim).*repmat(u(ii,:)',1,dim))/W_m_minus(ii);

end
MixtureUpdated.M=M;
MixtureUpdated.W=W_m_minus;
MixtureUpdated.Mu=Mu_minus;
MixtureUpdated.Sigma=Sigma_minus;