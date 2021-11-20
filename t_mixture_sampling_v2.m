function [X,tao]=t_mixture_sampling_v2(M,W_m,Mu,Sigma,df)
r=W_m(1,1);
rand_num=rand;
for ii=1:M % Component index of the Gaussian mixture
    if rand_num<=r
        tao=gamrnd(df(ii)/2,1/(df(ii)/2));% gamrnd(a,b)=gengamma(a,1/b);
        %X=mvnrnd(Mu(ii,:),Sigma(:,:,ii)/tao);
        X=mvnrnd(Mu(ii,:),Sigma(:,:,ii))/sqrt(tao);
        %cmp_ind=ii;
        break;
    else
        r=r+W_m(1,ii+1);
    end
end