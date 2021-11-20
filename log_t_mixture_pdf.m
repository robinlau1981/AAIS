function pdf =log_t_mixture_pdf(X,M,W,Mu,Sigma,df);%Ref: Robust mixture modelling using the t distribution; ML estimation of the t distribution using EM and its extensions,ECM and ECME
n=size(X,1);

log_Cpt_pdf=zeros(n,M);
parfor i=1:M
    log_Cpt_pdf(:,i)=log_t_pdf(X,Mu(i,:),Sigma(:,:,i),df);
end
pdf=log(exp(log_Cpt_pdf)*W');