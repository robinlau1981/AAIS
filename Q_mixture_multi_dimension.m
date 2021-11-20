function sum2=Q_mixture_multi_dimension(X,M,W_m,Mu,Sigma)
% sum=0;
% for i=1:M
%     sum=sum+W_m(1,i)*mvnpdf(X,Theta(1,i),Theta(2,i));
% end
sum2=0;
for i=1:M
    if sum(sum(Sigma(:,:,i)))>0
        sum2=sum2+W_m(1,i)*mvnpdf(X,Mu(i,:),Sigma(:,:,i));
    end
end