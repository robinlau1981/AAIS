function pdf =Gaussian_Mixture_pdf(X,M,W,Mu,Sigma);
sum2=0;
for i=1:M
    if cond(Sigma(:,:,i))>1e10
       %Sigma(:,:,i)=diag(ones(1,length(X)))*1e-2;
       Sigma(:,:,i)=Sigma(:,:,i)+.1*trace(Sigma(:,:,i))/length(X);
    end
    sum2=sum2+W(1,i)*mvnpdf(X,Mu(i,:),Sigma(:,:,i));
end
pdf=sum2;