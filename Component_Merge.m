% Component Merge Ref: Estimation of the number of components in a mixture
% model using stepwise split-and-merge EM Algorithm

clear;
load Res_SAIS_v1_02;
dim=length(Mu(1,:));
for k=1:M
    z(:,k)=mvnpdf(proposed_new_theta,Mu(k,:),Sigma(:,:,k));
end
for k=1:M
    z_mean(k)=mean(z(:,k));
end
cor=zeros(M,M);
% for i=1:M
%     for j=1:M
%         if i==j
%             cor(i,j)=1;
%         else
%             cor(i,j)=(z(:,i)-z_mean(i)*ones(N,1))'*(z(:,j)-z_mean(j)*ones(N,1))/(norm(z(:,i)-z_mean(i)*ones(N,1))*norm(z(:,j)-z_mean(j)*ones(N,1)));
%         end
%     end
% end
for i=1:M
    for j=i+1:M
        cor(i,j)=(z(:,i)-z_mean(i)*ones(N,1))'*(z(:,j)-z_mean(j)*ones(N,1))/(norm(z(:,i)-z_mean(i)*ones(N,1))*norm(z(:,j)-z_mean(j)*ones(N,1)));
    end
end
cor
cor_thr=.05; % Threshold
m_ind=ones(1,M);
for i=1:M
    if m_ind(i)==1
        for j=i+1:M
            if m_ind(j)==1
                if abs(cor(i,j))>cor_thr
                    m_ind(j)=0;                    
                    Mu(i,:)=W_m([i j])*Mu([i j],:)/sum(W_m([i j]));
                    Sigma(:,:,i)=(W_m(i)*Sigma(:,:,i)+W_m(j)*Sigma(:,:,j))/sum(W_m([i j]));
                    W_m(i)=sum(W_m([i j]));
                end
            end
        end
    end
end
M_remain_ind=find(m_ind==1);
M_remain=length(M_remain_ind);
Mu_new=zeros(M_remain,dim);
Sigma_new=zeros(dim,dim,M_remain);
W_m_new=zeros(1,M_remain);
for k=1:M_remain
    Mu_new(k,:)=Mu(M_remain_ind(k),:);
    Sigma_new(:,:,k)=Sigma(:,:,M_remain_ind(k));
    W_m_new(k)=W_m(M_remain_ind(k));
end
M=M_remain
Mu=Mu_new
Sigma=Sigma_new
W_m=W_m_new