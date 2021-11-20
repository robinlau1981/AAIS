function [Proposal,X,Sign]=Merge(Proposal,X,Cor_thr,df,dim)

%-------------Components Merging------------%
        Mixture_temp=Proposal;
        z=repmat(Mixture_temp.W,X.N,1).*X.Resp./repmat(X.Proposal,1,Proposal.M);
        %       z=zeros(N,M2);
        %         for k=1:Mixture_temp.M
        %             z(:,k)=W_m(k)*rp(:,k)./proposal; %z(:,k)=q.*mvnpdf(X,Mu2(k,:),Sigma2(:,:,k));
        %         end
        z_mean=X.NormalizedWeight'*z;
        cor=zeros(Mixture_temp.M,Mixture_temp.M);
        for r=1:Mixture_temp.M
            for c=r+1:Mixture_temp.M
                cor(r,c)=sum((z(:,r)-z_mean(r)).*(z(:,c)-z_mean(c)).*X.NormalizedWeight)/...
                    (sqrt(sum((z(:,r)-z_mean(r)).^2.*X.NormalizedWeight))*sqrt(sum((z(:,c)-z_mean(c)).^2.*X.NormalizedWeight)));
            end
        end
        
        m_ind=ones(1,Mixture_temp.M);
        for r=1:Mixture_temp.M
            if m_ind(r)==1
                for c=r+1:Mixture_temp.M
                    if m_ind(c)==1
                        if cor(r,c)>Cor_thr
                            m_ind(c)=0;
                            
                             ind_pts_c=find(X.Root==c);
%                             ind_pts_r=find(X.Root==r);
%                             X3.Values=[X.Values(ind_pts_c,:);X.Values(ind_pts_r,:)];
%                             X3.NormalizedWeight=X.NormalizedWeight([ind_pts_c;ind_pts_r]);
%                             X3.NormalizedWeight=X3.NormalizedWeight/sum(X3.NormalizedWeight);
%                             X3.N=length([ind_pts_c;ind_pts_r]);
%                             Mixture_temp.Mu(r,:)=X3.NormalizedWeight'*X3.Values;
%                             Mixture_temp.Sigma(:,:,r)=((X3.Values-repmat(Mixture_temp.Mu(r,:),X3.N,1)).*repmat(X3.NormalizedWeight,1,dim))'*(X3.Values-repmat(Mixture_temp.Mu(r,:),X3.N,1))*(df-2)/df;
%                             Mixture_temp.Sigma(:,:,r)=(Mixture_temp.Sigma(:,:,r)'+Mixture_temp.Sigma(:,:,r))/2;

                            Mixture_temp.Mu(r,:)=(Mixture_temp.Mu(r,:)*Mixture_temp.W(r)+Mixture_temp.Mu(c,:)*Mixture_temp.W(c))/(Mixture_temp.W(r)+Mixture_temp.W(c));
                            Mixture_temp.Sigma(:,:,r)=(Mixture_temp.Sigma(:,:,r)*Mixture_temp.W(r)+Mixture_temp.Sigma(:,:,c)*Mixture_temp.W(c))/(Mixture_temp.W(r)+Mixture_temp.W(c));
                            
                            if ~isposdef(Mixture_temp.Sigma(:,:,r))
                                if dim==7
                                    Mixture_temp.Sigma(:,:,r)=1e6*diag([1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6]);
                                elseif dim==12
                                    Mixture_temp.Sigma(:,:,r)=1e6*diag([1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6]);
                                elseif dim==17
                                    Mixture_temp.Sigma(:,:,r)=1e6*diag([1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6 1e-10 1e-6 1e-6 1e-6 1e-6]);
                                end
                            end
                            
                            
                            Mixture_temp.W(r)=sum(Mixture_temp.W([r c]));
                            X.Resp(:,r)=exp(log_t_pdf(X.Values,Mixture_temp.Mu(r,:),Mixture_temp.Sigma(:,:,r),df));
                            X.Root(ind_pts_c)=r*ones(length(ind_pts_c),1);
                        end
                    end
                end
            end
        end
        Sign=0;
        index_remain=find(m_ind==1);
        M_temp=length(index_remain);
        if M_temp<Proposal.M
            Proposal.M=length(index_remain);
            Proposal.Mu=Mixture_temp.Mu(index_remain,:);
            Proposal.Sigma=Mixture_temp.Sigma(:,:,index_remain);
            Proposal.W=Mixture_temp.W(index_remain);
            X.Resp=X.Resp(:,index_remain);
            X.Proposal=X.Resp*Proposal.W';
            X.logProposal=log(X.Proposal);
            X.logWeight=X.logAnnealTarget-X.logProposal;
            Root_temp=X.Root;
            for i=1:M_temp
                ind_temp=find(Root_temp==index_remain(i));
                X.Root(ind_temp)=i*ones(length(ind_temp),1);
            end
            logWeight_scaled=X.logWeight-max(X.logWeight)+10; % Scaling
            weight_temp=exp(logWeight_scaled);
            X.NormalizedWeight=weight_temp/sum(weight_temp);
            Sign=1;
            
        end