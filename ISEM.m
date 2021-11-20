function [X,Proposal] =ISEM(X,Proposal,Proposal0,df,Target_function,gama,j) % One iteration of EM with IS

    X=t_mix_sample(Proposal,X,df);
    index_survival=[];
    for i=1:Proposal.M
        if length(find(X.Root==i))>0%~(isempty(X.Root==i))
            index_survival=[index_survival i];
        end
    end
    if length(index_survival)<Proposal.M
        Proposal.M=length(index_survival);
        Proposal.Mu=Proposal.Mu(index_survival,:);
        Proposal.Sigma=Proposal.Sigma(:,:,index_survival);
        Proposal.W=Proposal.W(index_survival);
        Proposal.W=Proposal.W/sum(Proposal.W);
        
        Root_temp=X.Root;
        for i=1:length(index_survival)
            temp=find(Root_temp==index_survival(i));
            X.Root(temp)=i*ones(length(temp),1);
        end
    end
    
    X.Resp=zeros(X.N,Proposal.M); % responsibilily of each component with regard to each particle
    for i=1:Proposal.M
        X.Resp(:,i)=exp(log_t_pdf(X.Values,Proposal.Mu(i,:),Proposal.Sigma(:,:,i),df));
    end
    X.Proposal=X.Resp*Proposal.W';
    X.logProposal=log(X.Proposal);
    
    if j>1
        rp_temp=zeros(X.N,Proposal0.M);
        for i=1:Proposal0.M
            rp_temp(:,i)=exp(log_t_pdf(X.Values,Proposal0.Mu(i,:),Proposal0.Sigma(:,:,i),df));
        end
        X.logProposal0=log(rp_temp*Proposal0.W');
    else % j==1
        X.logProposal0=X.logProposal;
    end
    
    X.Target=feval(Target_function,X.Values);%Target_pdf_rastrigin(X.Values);%feval(Target_function,X.Values,data);
    X.logTarget=log(X.Target);
    
    X.logAnnealTarget=X.logTarget*gama(j)+X.logProposal0*(1-gama(j));
       
    X.logWeight=X.logAnnealTarget-X.logProposal;
    logWeight_scaled=X.logWeight-max(X.logWeight)+10; % Scaling
    weight_temp=exp(logWeight_scaled);
    X.NormalizedWeight=weight_temp/sum(weight_temp);
    
%     ESS(j)=1/sum(X.NormalizedWeight.^2)/X.N;
%     fprintf('ESS before EM: %i \r', ESS(j));
%     str1=num2str(j);
%     str2='th iter--Before EM';
%     strf = [str1 str2 str0];
%     
%     Plot_GM(X.Values(:,1:2),Proposal.M,Proposal.W,Proposal.Mu(:,1:2)',df/(df-2)*Proposal.Sigma(1:2,1:2,:));text(text_Pos_X,text_Pos_Y,strf,'BackgroundColor',[.7 .9 .7]);axis(axis_lmt);grid on;
%     aviobj = addframe(aviobj, getframe);
    
    Proposal=t_mix_update_v2(X,Proposal,df);
    Proposal.W=Proposal.W/sum(Proposal.W);
    X=t_mix_sample(Proposal,X,df);
    
    index_survival=[];
    for i=1:Proposal.M
        if length(find(X.Root==i))>0%~(isempty(X.Root==i))
            index_survival=[index_survival i];
        end
    end
    if length(index_survival)<Proposal.M
        Proposal.M=length(index_survival);
        Proposal.Mu=Proposal.Mu(index_survival,:);
        Proposal.Sigma=Proposal.Sigma(:,:,index_survival);
        Proposal.W=Proposal.W(index_survival);
        Proposal.W=Proposal.W/sum(Proposal.W);
        
        Root_temp=X.Root;
        for i=1:length(index_survival)
            temp=find(Root_temp==index_survival(i));
            X.Root(temp)=i*ones(length(temp),1);
        end
    end
    
    X.Resp=zeros(X.N,Proposal.M); % responsibilily of each component with regard to each particle
    for i=1:Proposal.M
        X.Resp(:,i)=exp(log_t_pdf(X.Values,Proposal.Mu(i,:),Proposal.Sigma(:,:,i),df));
    end
    X.Proposal=X.Resp*Proposal.W';
    X.logProposal=log(X.Proposal);
    
    rp_temp=zeros(X.N,Proposal0.M);
    for i=1:Proposal0.M
        rp_temp(:,i)=exp(log_t_pdf(X.Values,Proposal0.Mu(i,:),Proposal0.Sigma(:,:,i),df));
    end
    X.logProposal0=log(rp_temp*Proposal0.W');
    X.Target=Target_pdf_rastrigin(X.Values);%X.Target=feval(Target_function,X.Values,data);
    X.logTarget=log(X.Target);
    %X.logAnnealTarget=X.logPrior+X.logLike*gama(j);
    X.logAnnealTarget=X.logTarget*gama(j)+X.logProposal0*(1-gama(j));
    
    X.logWeight=X.logAnnealTarget-X.logProposal;
    logWeight_scaled=X.logWeight-max(X.logWeight)+10; % Scaling
    weight_temp=exp(logWeight_scaled);
    X.NormalizedWeight=weight_temp/sum(weight_temp);
    
    %ESS(j)=1/sum(X.NormalizedWeight.^2)/X.N;