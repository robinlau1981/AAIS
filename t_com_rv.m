clear;
close all;

load -ASCII hd73526_v2;
data.t=hd73526_v2(:,1);
data.V=hd73526_v2(:,2);
data.errors=hd73526_v2(:,3);

t=[min(data.t)-100:10:max(data.t)+100];
load Res_SAIS_v12_v7_2d_0p_v2_data2_02;
rv_0p=rv_final;
load Res_SAIS_v12_v7_7d_v2_data2;
N=1e5;
dim=7;
Y=zeros(N,dim);
Yt=Y;
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
end
q=Target_pdf_v6(Y,data)./t_mixture_pdf(Y,M,W_m,Mu,Sigma,df);
MarLik=mean(q);
MarLik
s0_square=mean((q-MarLik).^2);
ci=3*sqrt(s0_square)/sqrt(N)
q=q/sum(q);
ESS_final=1/sum(q.^2)

str1 = 'title(';
str2 = ['Particle Size: N=' num2str(N)];
strf = [str1 'str2' ')'];
eval(strf);

%figure,plot(gama,ESS/N);xlabel('\phi');ylabel('ESS/N');

eval(strf);

rv_1p=q'*rv_est;
% load Res_SAIS_v12_v7_12d_v2_data2_02;
% rv_2p=rv_final;
% load Res_SAIS_v12_v7_12d_v2_data2;
% rv_2p=rv_final;
load Res_v12_12d_trans_para_v4_data2_02;
dim=12;
Yt=zeros(N,dim);
for i=1:N    
    Yt(i,:)=Para_transform_12d(Y(i,:));
    rv_est(i,:)=(model_v4(Yt(i,:),t))';
end

q=Target_pdf_v7(Y,data)./t_mixture_pdf_v2(Y,M,W_m,Mu,Sigma,df);
q=q/sum(q);
rv_2p=q'*rv_est;

figure,
errorbar(data.t,data.V,data.errors);hold on;
plot(t,rv_0p,'r',t,rv_1p,'-.',t,rv_2p,'-.');xlabel('days');ylabel('RV');
hold off;

save Res_rv_comparison_data2 rv_0p rv_1p rv_2p;