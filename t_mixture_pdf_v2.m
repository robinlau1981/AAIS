function pdf =t_mixture_pdf_v2(X,M,W,Mu,Sigma,df) %Ref: Robust mixture modelling using the t distribution; ML estimation of the t distribution using EM and its extensions,ECM and ECME
[n,p]=size(X);
pdf=zeros(n,1);
%sum2=pdf;
% for i=1:M
%     sum2=sum2+W(1,i).*t_pdf(X,Mu(i,:),Sigma(:,:,i),df);
% end
Cpt_pdf=zeros(n,M);
for i=1:M
    Cpt_pdf(:,i)=t_pdf(X,Mu(i,:),Sigma(:,:,i),df(i));
end
pdf=Cpt_pdf*W';