function X=t_mix_sample(Proposal,X,df)
M=Proposal.M;
W_m=Proposal.W;
Mu=Proposal.Mu;
Sigma=Proposal.Sigma;
dim=size(Mu,2);
pts=zeros(X.N,dim);
root=zeros(X.N,1);
%parfor i=1:X.N  % Particle index
for i=1:X.N  % Particle index
    [pts(i,:),root(i)]=t_mixture_sampling(M,W_m,Mu,Sigma,df);
end
X.Values=pts;
X.Root=root;