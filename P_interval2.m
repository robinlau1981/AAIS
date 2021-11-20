clear
D=input('which data set \n');
M=input('which model \n');

if D==1
    load -ASCII hd73526;
    data.t=hd73526(:,1);
    data.V=hd73526(:,2);
    data.errors=hd73526(:,3);
else
    load -ASCII hd73526_v2;
    data.t=hd73526_v2(:,1);
    data.V=hd73526_v2(:,2);
    data.errors=hd73526_v2(:,3);
end

matlabpool open;
N=5e5;
if M==1    
    if D==1
        load Res_7d_trans_para_data1;
    else
        load Res_7d_trans_para_data2;
    end
    
    
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
else
    % For 2p model
    if D==1
        load Res_12d_trans_para_data1;
    else
        load Res_v12_12d_trans_para_v4_data2_02;
    end
    dim=12;
    Y=zeros(N,dim);
    
    parfor i=1:N
        Y(i,:)=t_mixture_sampling_v2(M,W_m,Mu,Sigma,df);
    end
    P_1p=exp(Y(:,3));
    P_2p=exp(Y(:,8));
    q=Target_pdf_v7(Y,data)./t_mixture_pdf_v2(Y,M,W_m,Mu,Sigma,df);
    
    q=q/sum(q);
    ESS=1/sum(q.^2)
    % For data1
    [P_sort,ind]=sort(P_1p);
    q_sort=q(ind);
    
    N1=length(find(P_sort<300));
    ratio1=sum(q_sort(1:N1))/sum(q_sort(N1+1:N))
    
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
    
    [P_sort,ind]=sort(P_2p);
    q_sort=q(ind);
    
    N1=length(find(P_sort<300));
    ratio2=sum(q_sort(1:N1))/sum(q_sort(N1+1:N))
    
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
end
matlabpool close;