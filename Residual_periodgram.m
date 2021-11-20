load -ASCII hd73526_v2;
data.t=hd73526_v2(:,1);
data.V=hd73526_v2(:,2);
data.errors=hd73526_v2(:,3);

t=data.t;
MarLik=0;
q=zeros(N,1);
rv_est=zeros(N,length(t));
%rv_final=q*rv_est;
for i=1:N
    r=W_m(1,1);
    rand_num=rand;
    for ii=1:M
        if rand_num<=r
            tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
            Y(i,:)=mvnrnd(Mu(ii,:),Sigma(:,:,ii)/tao);
            break;
        else
            r=r+W_m(1,ii+1);
        end
    end

    Yt(i,:)=Y(i,:);
    Yt(i,2)=exp(Y(i,2));
    Yt(i,3)=exp(Y(i,3));
    Yt(i,4)=sqrt(Y(i,4)^2+Y(i,5)^2);
    Yt(i,5)=mod(atan(Y(i,5)/Y(i,4))+pi*(Y(i,4)<0),2*pi);
    Yt(i,6)=mod(Y(i,6)-Yt(i,5),2*pi);
    Yt(i,7)=exp(Y(i,7));
    for j=1:length(t)
        rv_est(i,j)=model_v2(Yt(i,1:6),t(j));
    end
    Resi_rv=data.V-rv_est(i,:)';
    V=Resi_rv;
    omega=1:.1:1000;
    Scargle_omega=periodogram_scargle(t,V,omega);
    figure,plot(omega,Scargle_omega);xlabel('period candidates (days)');
end
proposal=t_mixture_pdf(Y,M,W_m,Mu,Sigma,df);%Gaussian_Mixture_pdf(X(i,:),M,W_m,Mu,Sigma);
q=Target_pdf_v6(Y,data)./proposal;

MarLik=mean(q);
MarLik
s0_square=mean((q-MarLik).^2);
ci=3*sqrt(s0_square)/sqrt(N)
q=q/sum(q);

rv_final=q'*rv_est;
Resi_rv=data.V-rv_final';


V=Resi_rv;
omega=1:.1:1000;
Scargle_omega=periodogram_scargle(t,V,omega);
semilogx(omega,Scargle_omega);xlabel('period candidates (days)');title('Periodogram for 51 Peg data');