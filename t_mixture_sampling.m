function [X,cmp_ind]=t_mixture_sampling(M,W_m,Mu,Sigma,df)
r=W_m(1,1);
dim=length(Mu);
rand_num=rand;
for ii=1:M % Component index of the Gaussian mixture
    if rand_num<=r
        tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
        %X=mvnrnd(Mu(ii,:),Sigma(:,:,ii)/tao);
        Sigma_temp=Sigma(:,:,ii)/tao;
        if isposdef(Sigma_temp)
            X=(randnorm(1,Mu(ii,:)',[],Sigma_temp))';
        else
            eig_values=eig(Sigma_temp);
            delta=unifrnd(mean(eig_values),max(eig_values));
            Sigma_temp=Sigma_temp+1e-10*delta*eye(dim);
            X=(randnorm(1,Mu(ii,:)',[],Sigma_temp))';
        end
        cmp_ind=ii;
        break;
    else
        r=r+W_m(1,ii+1);
    end
end