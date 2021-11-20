function MixtureUpdated=t_mix_update_v4(X,Proposal,df)

M=Proposal.M;
W_m=Proposal.W;
Mu=Proposal.Mu;
Sigma=Proposal.Sigma;
q=X.NormalizedWeight;
rp=X.Resp;
proposal=X.Proposal;
XV=X.Values;

[N,dim]=size(X.Values);

W_m_minus=zeros(1,M);
Mu_minus=zeros(M,dim);
Sigma_minus=zeros(dim,dim,M);

vz=dim;
if dim==17
    Sigma_z=diag([1e-2 1e-4 1e-10 1e-6 1e-6 1e-6 1e-4 1e-10 1e-6 1e-6 1e-6 1e-4 1e-10 1e-6 1e-6 1e-6 1e-2]);
elseif dim==12
    Sigma_z=diag([1e-2 1e-4 1e-10 1e-6 1e-6 1e-6 1e-4 1e-10 1e-6 1e-6 1e-6 1e-2]);
elseif dim==7
    Sigma_z=diag([1e-2 1e-4 1e-10 1e-6 1e-6 1e-6 1e-2]);
else
    Sigma_z=eye(dim)*1e-6;
end
    
    parfor ii=1:M
        rou(:,ii)=W_m(ii)*rp(:,ii)./proposal;
        W_m_minus(ii)=rou(:,ii)'*q;
        u(:,ii)=(df+dim)*ones(N,1)./(df*ones(N,1)+sqdist(XV',Mu(ii,:)',inv_posdef(Sigma(:,:,ii))));%MahalanobisDist(X.Values,,Sigma(:,:,ii)));%sqdist(X.Values',Mu(ii,:)',inv_posdef(Sigma(:,:,ii)))');
        weight_temp=(rou(:,ii).*u(:,ii).*q)/((rou(:,ii).*u(:,ii))'*q);
        Mu_minus(ii,:)=weight_temp'*XV;
        Nz=(rou(:,ii).*u(:,ii))'*q*N;
        weight_temp2=rou(:,ii).*u(:,ii).*q*N; % weight_temp2=weight_temp*Nz
        matrix_temp=XV-repmat(Mu_minus(ii,:),N,1);
        %Sigma_minus(:,:,ii)=matrix_temp'*(matrix_temp.*repmat(weight_temp,1,dim));
        %Sigma_minus(:,:,ii)=matrix_temp'*(matrix_temp.*repmat(weight_temp2,1,dim))/Nz;
        Sigma_minus(:,:,ii)=matrix_temp'*(matrix_temp.*repmat(weight_temp2,1,dim))/(Nz+(vz+dim+1))+vz*Sigma_z/(Nz+(vz+dim+1));
        Sigma_minus(:,:,ii)=(Sigma_minus(:,:,ii)'+Sigma_minus(:,:,ii))/2;
       
    end

MixtureUpdated.M=M;
MixtureUpdated.W=W_m_minus/sum(W_m_minus);
MixtureUpdated.Mu=Mu_minus;
MixtureUpdated.Sigma=Sigma_minus;