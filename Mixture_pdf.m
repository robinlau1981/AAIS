function pdf =Mixture_pdf(X,M,W,Mu,Sigma);
sum2=0;
for i=1:M
    sum2=sum2+W(1,i)*mvnpdf(X,Mu(i,:),Sigma(:,:,i));
end
pdf=sum2;