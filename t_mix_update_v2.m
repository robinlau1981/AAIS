function MixtureUpdated=t_mix_update_v2(X,Proposal,df)

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
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% rou=zeros(M,N);
% u=zeros(M,N);

parfor ii=1:M
    rou(:,ii)=W_m(ii)*rp(:,ii)./proposal;
    W_m_minus(ii)=rou(:,ii)'*q;
    u(:,ii)=(df+dim)*ones(N,1)./(df*ones(N,1)+sqdist(XV',Mu(ii,:)',inv_posdef(Sigma(:,:,ii))));%MahalanobisDist(X.Values,,Sigma(:,:,ii)));%sqdist(X.Values',Mu(ii,:)',inv_posdef(Sigma(:,:,ii)))');
    weight_temp=(rou(:,ii).*u(:,ii).*q)/((rou(:,ii).*u(:,ii))'*q);
    Mu_minus(ii,:)=weight_temp'*XV;
    matrix_temp=XV-repmat(Mu_minus(ii,:),N,1);
    Sigma_minus(:,:,ii)=matrix_temp'*(matrix_temp.*repmat(weight_temp,1,dim));
    Sigma_minus(:,:,ii)=(Sigma_minus(:,:,ii)'+Sigma_minus(:,:,ii))/2;
    if ~isposdef(Sigma_minus(:,:,ii))
        if dim==7
            Sigma_minus(:,:,ii)=1e6*diag([1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6]);
        elseif dim==2
            Sigma_minus(:,:,ii)=1e6*diag([1e-6 1e-6]);
        elseif dim==12
            Sigma_minus(:,:,ii)=1e6*diag([1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6]);
        else
            Sigma_minus(:,:,ii)=1e6*diag(1e-6*ones(1,dim));
        end
        %X_temp=XV+(randnorm(N,zeros(dim,1),[],Sigma_temp))';
        %matrix_temp=X_temp-repmat(Mu_minus(ii,:),N,1);
        %matrix_temp'*(matrix_temp.*repmat(weight_temp,1,dim));
        Sigma_minus(:,:,ii)=(Sigma_minus(:,:,ii)'+Sigma_minus(:,:,ii))/2;
    end
    %eig_values=eig(Sigma_minus(:,:,ii));
    %delta=(mean(eig_values)+max(eig_values))/2;
    %Sigma_minus(:,:,ii)=Sigma_minus(:,:,ii)+1e-10*delta*eye(size(X,2));
    %Sigma_minus(:,:,ii)=Sigma_minus(:,:,ii)+1e-10*eye(size(X.Values,2));
end
MixtureUpdated.M=M;
MixtureUpdated.W=W_m_minus;
MixtureUpdated.Mu=Mu_minus;
MixtureUpdated.Sigma=Sigma_minus;