function X=t_mixture_sampling_v3(M,W_m,Mu,Sigma,df)
r=W_m(1,1);
rand_num=rand;
dim=size(Mu,2);
for ii=1:M % Component index of the Gaussian mixture
    if rand_num<=r
        tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
%         if rcond(Sigma(:,:,ii)/tao)<1e-5
%             pause;
%         end
        X=Mu(ii,:)+randn(1,dim)*chol(Sigma(:,:,ii)/tao);
        %Sigma_s=Sigma(:,:,ii);
        break;
    else
        r=r+W_m(1,ii+1);
    end
end

% for ii=1:2e5
%     gr(ii)=gamrnd(df/2,1/(df/2));
% end
% figure,plot(gr);