clear
load -ASCII hd73526;
data.t=hd73526(:,1);
data.V=hd73526(:,2);
data.errors=hd73526(:,3);

% load -ASCII hd73526_v2;
% data.t=hd73526_v2(:,1);
% data.V=hd73526_v2(:,2);
% data.errors=hd73526_v2(:,3);

load Res_7d_trans_para_data1;
%load Res_7d_trans_para_data2;
matlabpool open;
N=1e5;
dim=7;
Y=zeros(N,dim);
Yt=zeros(N,dim);
P=zeros(N,1);
parfor i=1:N
    Y(i,:)=t_mixture_sampling_v2(M,W_m,Mu,Sigma,df);    
end
P=exp(Y(:,3));
q=Target_pdf_v6(Y,data)./t_mixture_pdf_v2(Y,M,W_m,Mu,Sigma,df);

q=q/sum(q); 
Ytminus=zeros(N,dim);
% P=zeros(N,1);
% 
% for i = 1 : N
%     P(i)=unif_resample(Yt(:,3),q);    
% end


[P_sort,ind]=sort(P);
q_sort=q(ind);

N1=length(find(P_sort<300));
ratio=sum(q_sort(1:N1))/sum(q_sort(N1+1:N))

q_1=q_sort(1:N1)/sum(q_sort(1:N1));
q_2=q_sort(N1+1:N)/sum(q_sort(N1+1:N));
P_1=P_sort(1:N1);
P_2=P_sort(N1+1:N);

P_1_mean=sum(q_1.*P_1)
lb=P_1(length(find(cumsum(q_1)<.025)))
ub=P_1(length(find(cumsum(q_1)<.975)))

P_2_mean=sum(q_2.*P_2)
lb=P_2(length(find(cumsum(q_2)<.025)))
ub=P_2(length(find(cumsum(q_2)<.975)))


% i_1=(P<300);
% i_2=(P>300);
% q_1=q(find(i_1==1));
% q_2=q(find(i_2==1));
% 
% ratio=length(q_1)/length(q_2)
% 
% L1=length(q_1)
% L2=length(q_2)
% 
% q_1=q_1/sum(q_1);
% P_1_mean=sum(q_1.*P(i_1))
% P1=sort(P(i_1));
% lb=P1(ceil(0.025*L1))
% ub=P1(floor(0.975*L1))
% 
% %s_1=sum(q_1.*((P(i_1)-P_1_mean).^2));
% %ci_1=3*sqrt(s_1)%/sqrt(length(find(i_1==1)))
% 
% q_2=q_2/sum(q_2);
% P_2_mean=sum(q_2.*P(i_2))
% P2=sort(P(i_2));
% lb2=P2(ceil(0.025*L2))
% ub2=P2(ceil(0.975*L2))


%s_2=sum(q_2.*((P(i_2)-P_2_mean).^2));
%ci_2=3*sqrt(s_2)%/sqrt(length(find(i_2==1)))

matlabpool close;

%load Res_12d_trans_para_data1;
% load Res_v12_12d_trans_para_v4_data2_02;
% matlabpool open;
% N=5e5;
% dim=12;
% Y=zeros(N,12);
% P_1p=zeros(N,1);
% P_2p=zeros(N,1);
% parfor i=1:N
%     Y(i,:)=t_mixture_sampling_v2(M,W_m,Mu,Sigma,df);
%     %Yt(i,:)=Para_transform_12d(Y(i,:));
% end
% P_1p=exp(Y(:,3));
% P_2p=exp(Y(:,8));
% q=Target_pdf_v7(Y,data)./t_mixture_pdf_v2(Y,M,W_m,Mu,Sigma,df);
% 
% q=q/sum(q); 
% % Ytminus=zeros(N,dim);
% % 
% % for i = 1 : N
% %     P_1p(i)=unif_resample(Yt(:,3),q);
% %     P_2p(i)=unif_resample(Yt(:,3+5),q);
% % end
% 
% [P_sort,ind]=sort(P_1p);
% q_sort=q(ind);
% 
% N1=length(find(P_sort<300));
% ratio1=sum(q_sort(1:N1))/sum(q_sort(N1+1:N))
% 
% q_1=q_sort(1:N1)/sum(q_sort(1:N1));
% q_2=q_sort(N1+1:N)/sum(q_sort(N1+1:N));
% P_1=P_sort(1:N1);
% P_2=P_sort(N1+1:N);
% 
% P_1_mean=sum(q_1.*P_1)
% lb=P_1(length(find(cumsum(q_1)<.025)))
% ub=P_1(length(find(cumsum(q_1)<.975)))
% 
% P_2_mean=sum(q_2.*P_2)
% lb=P_2(length(find(cumsum(q_2)<.025)))
% ub=P_2(length(find(cumsum(q_2)<.975)))
% 
% [P_sort,ind]=sort(P_2p);
% q_sort=q(ind);
% 
% N1=length(find(P_sort<300));
% ratio2=sum(q_sort(1:N1))/sum(q_sort(N1+1:N))
% 
% q_1=q_sort(1:N1)/sum(q_sort(1:N1));
% q_2=q_sort(N1+1:N)/sum(q_sort(N1+1:N));
% P_1=P_sort(1:N1);
% P_2=P_sort(N1+1:N);
% 
% P_1_mean=sum(q_1.*P_1)
% lb=P_1(length(find(cumsum(q_1)<.025)))
% ub=P_1(length(find(cumsum(q_1)<.975)))
% 
% P_2_mean=sum(q_2.*P_2)
% lb=P_2(length(find(cumsum(q_2)<.025)))
% ub=P_2(length(find(cumsum(q_2)<.975)))

% i_1p_1=(P_1p<300);
% i_1p_2=(P_1p>300);
% q_1p_1=q(i_1p_1);
% q_1p_2=q(i_1p_2);
% 
% i_2p_1=(P_2p<300);
% i_2p_2=(P_2p>300);
% q_2p_1=q(i_2p_1);
% q_2p_2=q(i_2p_2);
% 
% q_1p_1=q_1p_1/sum(q_1p_1);
% P_1p_1_mean=sum(q_1p_1.*P_1p(i_1p_1))
% s_1p_1=sum(q_1p_1.*((P_1p(i_1p_1)-P_1p_1_mean).^2));
% ci_1p_1=3*sqrt(s_1p_1)%/sqrt(length(find(i_1p_1==1)))
% 
% q_1p_2=q_1p_2/sum(q_1p_2);
% P_1p_2_mean=sum(q_1p_2.*P_1p(i_1p_2))
% s_1p_2=sum(q_1p_2.*((P_1p(i_1p_2)-P_1p_2_mean).^2));
% ci_1p_2=3*sqrt(s_1p_2)%/sqrt(length(find(i_1p_2==1)))
% 
% q_2p_1=q_2p_1/sum(q_2p_1);
% P_2p_1_mean=sum(q_2p_1.*P_2p(i_2p_1))
% s_2p_1=sum(q_2p_1.*((P_2p(i_2p_1)-P_2p_1_mean).^2));
% ci_2p_1=3*sqrt(s_2p_1)%/sqrt(length(find(i_2p_1==1)))
% 
% q_2p_2=q_2p_2/sum(q_2p_2);
% P_2p_2_mean=sum(q_2p_2.*P_2p(i_2p_2))
% s_2p_2=sum(q_2p_2.*((P_2p(i_2p_2)-P_2p_2_mean).^2));
% ci_2p_2=3*sqrt(s_2p_2)%/sqrt(length(find(i_2p_2==1)))

matlabpool close;