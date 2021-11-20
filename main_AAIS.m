% Code for implementing the adaptive annealed importance sampling algorithm
% published in the paper "Liu, B., Adaptive Annealed Importance Sampling for Multimodal Posterior Exploration and Model Selection with Application to Extrasolar Planet Detection, The Astrophysical Journal Supplement Series,
% vol. 213, no.14, pp.1-16, 2014. doi:10.1088/0067-0049/213/1/14."
clear;
close all;
addpath(genpath('D:\Code\Matlab\work\Toolbox_plus\lightspeed')); % add folder address where Tom Minka's lightspeed toolbox is saved
Target_function='work_3p_transformed'; % Target function associated with the posterior, i.e. the target distribution to be sampled; 
% 'work_3p_transformed' corresponds to the 3-planet model in the paper

% %---------------------------Simulating Annealing Importance Sampling-------------------------------%

switch Target_function
    case 'work_3p_transformed'
        load -ASCII hd73526_v2;
        data.t=hd73526_v2(:,1);
        data.V=hd73526_v2(:,2);
        data.errors=hd73526_v2(:,3);
        dim=17;
        
        N=4e5; X2.N=1e4; % N: particle size of importance sampling; X2.N: sample size for adding a new mixure component
        
        min_C=-200; max_C=200; min_K=1; max_K=200; min_P=1;max_P=2000;min_s=1;max_s=200; %define value spaces of the model parameters
        
        Proposal.M=50;  % initial number of mixture components in the Proposal distribution
        Proposal.W=ones(1,Proposal.M)*1/Proposal.M; %  initial weights of mixture components in the Proposal distribution
        Proposal.Mu=zeros(Proposal.M,dim); %  initial centers of mixture components in the Proposal distribution
        Proposal.Mu(:,1)=unifrnd(min_C,max_C,Proposal.M,1);
        Proposal.Mu(:,2)=unifrnd(log(min_K),log(max_K),Proposal.M,1); % logK
        Proposal.Mu(:,3)=unifrnd(log(min_P),log(max_P),Proposal.M,1); % logP
        e=unifrnd(0,1,Proposal.M,1);
        w=unifrnd(0,2*pi,Proposal.M,1);
        M0=unifrnd(0,2*pi,Proposal.M,1);
        Proposal.Mu(:,4)=e.*cos(w);
        Proposal.Mu(:,5)=e.*sin(w); % y
        Proposal.Mu(:,6)=mod(w+M0,2*pi);
        Proposal.Mu(:,7)=unifrnd(log(min_K),log(max_K),Proposal.M,1); % logK
        Proposal.Mu(:,8)=unifrnd(log(min_P),log(max_P),Proposal.M,1); % logP
        e=unifrnd(0,1,Proposal.M,1);
        w=unifrnd(0,2*pi,Proposal.M,1);
        M0=unifrnd(0,2*pi,Proposal.M,1);
        Proposal.Mu(:,9)=e.*cos(w);
        Proposal.Mu(:,10)=e.*sin(w); % y
        Proposal.Mu(:,11)=mod(w+M0,2*pi);
        Proposal.Mu(:,12)=unifrnd(log(min_K),log(max_K),Proposal.M,1); % logK
        Proposal.Mu(:,13)=unifrnd(log(min_P),log(max_P),Proposal.M,1); % logP
        e=unifrnd(0,1,Proposal.M,1);
        w=unifrnd(0,2*pi,Proposal.M,1);
        M0=unifrnd(0,2*pi,Proposal.M,1);
        Proposal.Mu(:,14)=e.*cos(w);
        Proposal.Mu(:,15)=e.*sin(w); % y
        Proposal.Mu(:,16)=mod(w+M0,2*pi);
        Proposal.Mu(:,17)=unifrnd(log(min_s),log(max_s),Proposal.M,1); % logs
        Proposal.W=1/Proposal.M*ones(1,Proposal.M);
        df=5; % degree of freedom of the student's t distribution
        Proposal.Sigma=zeros(dim,dim,Proposal.M);
        fprintf('hd73526_v2 \r');
        gama=[.01 .1:.1:1]; % annealing schedule
end

temp=zeros(dim,dim);
temp(1,1)=1000;
temp(2,2)=10;
temp(3,3)=10;
temp(4,4)=.1;
temp(5,5)=.1;
temp(6,6)=10;
temp(7,7)=10;
temp(8,8)=10;
temp(9,9)=.1;
temp(10,10)=.1;
temp(11,11)=10;
temp(12,12)=10;
temp(13,13)=10;
temp(14,14)=.1;
temp(15,15)=.1;
temp(16,16)=10;
temp(17,17)=10;

Sigma_initial=temp/10;
for i=1:Proposal.M
    Proposal.Sigma(:,:,i)=temp;
end

Proposal0=Proposal;
ESS_Arr=[];

Cor_thr=.6;
ESS_thr=.6;
str0='-----Made by Bin Liu-----';
plot_id=0;
w_thr=1e-5; % If the weights of a mixture component is smaller than this value, delete it
for j=1:length(gama)% while gama(end)<=1
    X.N=N;
    X.Values=zeros(X.N,dim);
    X.Root=zeros(X.N,1);
    % One iteration of IS with EM
    if j==1
        %%  sampling from uniform distr
        X.Values(:,1)=unifrnd(min_C,max_C,X.N,1);
        X.Values(:,2)=unifrnd(log(min_K),log(max_K),X.N,1);
        X.Values(:,3)=unifrnd(log(min_P),log(max_P),X.N,1);
        e=unifrnd(0,1,X.N,1);
        w=unifrnd(0,2*pi,X.N,1);
        M0=unifrnd(0,2*pi,X.N,1);
        X.Values(:,4)=e.*cos(w);
        X.Values(:,5)=e.*sin(w);
        X.Values(:,6)=mod(w+M0,2*pi);
        X.Values(:,7)=unifrnd(log(min_K),log(max_K),X.N,1);
        X.Values(:,8)=unifrnd(log(min_P),log(max_P),X.N,1);
        e=unifrnd(0,1,X.N,1);
        w=unifrnd(0,2*pi,X.N,1);
        M0=unifrnd(0,2*pi,X.N,1);
        X.Values(:,9)=e.*cos(w);
        X.Values(:,10)=e.*sin(w);
        X.Values(:,11)=mod(w+M0,2*pi);
        X.Values(:,12)=unifrnd(log(min_K),log(max_K),X.N,1);
        X.Values(:,13)=unifrnd(log(min_P),log(max_P),X.N,1);
        e=unifrnd(0,1,X.N,1);
        w=unifrnd(0,2*pi,X.N,1);
        M0=unifrnd(0,2*pi,X.N,1);
        X.Values(:,14)=e.*cos(w);
        X.Values(:,15)=e.*sin(w);
        X.Values(:,16)=mod(w+M0,2*pi);
        X.Values(:,17)=unifrnd(log(min_s),log(max_s),X.N,1);
        %%
        
        X.Resp=zeros(X.N,Proposal.M); % responsibilily of each component with regard to each particle
        for i=1:Proposal.M
            X.Resp(:,i)=exp(log_t_pdf(X.Values,Proposal.Mu(i,:),Proposal.Sigma(:,:,i),df));
        end
        X.Proposal=X.Resp*Proposal.W';
        X.logProposal=log(X.Proposal);
        X.logProposal0=X.logProposal;
        
        
        [X.logPrior,X.logLike]=feval(Target_function,X.Values,data);
        X.logTarget=X.logPrior+X.logLike;
        X.logAnnealTarget=X.logTarget*gama(j)+X.logProposal0*(1-gama(j));
        
        X.logWeight=X.logAnnealTarget-X.logProposal;
        logWeight_scaled=X.logWeight-max(X.logWeight)+10; % Scaling
        weight_temp=exp(logWeight_scaled);
        X.NormalizedWeight=weight_temp/sum(weight_temp);
        
        Proposal=t_mix_update_v2(X,Proposal,df);
        Proposal.W=Proposal.W/sum(Proposal.W);
        X=t_mix_sample(Proposal,X,df);
        
        
        X.Resp=zeros(X.N,Proposal.M); % responsibilily of each component with regard to each particle
        for i=1:Proposal.M
            X.Resp(:,i)=exp(log_t_pdf(X.Values,Proposal.Mu(i,:),Proposal.Sigma(:,:,i),df));
        end
        X.Proposal=X.Resp*Proposal.W';
        X.logProposal=log(X.Proposal);
        
        rp_tmp=zeros(X.N,Proposal0.M);
        for i=1:Proposal0.M
            rp_tmp(:,i)=exp(log_t_pdf(X.Values,Proposal0.Mu(i,:),Proposal0.Sigma(:,:,i),df));
        end
        X.logProposal0=log(rp_tmp*Proposal0.W');
        
        [X.logPrior,X.logLike]=feval(Target_function,X.Values,data);
        X.logTarget=X.logPrior+X.logLike;
        X.logAnnealTarget=X.logTarget*gama(j)+X.logProposal0*(1-gama(j));
        
        X.logWeight=X.logAnnealTarget-X.logProposal;
        logWeight_scaled=X.logWeight-max(X.logWeight)+10; % Scaling
        weight_temp=exp(logWeight_scaled);
        X.NormalizedWeight=weight_temp/sum(weight_temp);
    else
        [X,Proposal] =ISEM_astro(X,Proposal,Proposal0,df,Target_function,gama,j,data);
    end
    
    
    ESS(j)=1/sum(X.NormalizedWeight.^2)/X.N;
    fprintf('ESS after EM: %i \r', ESS(j));
    
    
    str1=num2str(j);
    str2='th iter--After EM';
    strf = [str1 str2 str0];
    if plot_id
        Plot_GM(X.Values(:,1:2),Proposal.M,Proposal.W,Proposal.Mu(:,1:2)',df/(df-2)*Proposal.Sigma(1:2,1:2,:));axis(axis_lmt);grid on;
        aviobj = addframe(aviobj, getframe);
    end
    [Proposal,X,Sign]=Merge(Proposal,X,Cor_thr,df,dim);
    if Sign==1
        ESS(j)=1/sum(X.NormalizedWeight.^2)/X.N;
        %--------------------------------
        str1=num2str(j);
        str2='th iter--After Merging--ESS=';
        str3=' %i \r';
        strf = [str1 str2 str3];
        fprintf(strf, ESS(j));
        
        str2='th iter--After Merging';
        strf = [str1 str2 str0];
        if plot_id
            Plot_GM(X.Values(:,1:2),Proposal.M,Proposal.W,Proposal.Mu(:,1:2)',df/(df-2)*Proposal.Sigma(1:2,1:2,:));axis(axis_lmt);grid on; %text(text_Pos_X,text_Pos_Y,strf,'BackgroundColor',[.7 .9 .7]);
            aviobj = addframe(aviobj, getframe);
        end
    end
    
    
    %----------add component---------%
    j2=1;
    ESS_add=ESS(j);
    while ESS(j)<ESS_thr
        rp_temp=zeros(X.N,2);
        [q_max,ind_pt]=max(X.NormalizedWeight);
        ind_pt=ind_pt(1); % index of maximum weight particle
        Proposal_temp.Mu(1,:)=X.Values(ind_pt,:);
        Proposal_temp.Mu(2,:)=-Proposal_temp.Mu(1,:); %
        Proposal_temp.Sigma(:,:,1)=Sigma_initial;
        Proposal_temp.Sigma(:,:,2)=Sigma_initial;
        Proposal_temp.M=2;
        Proposal_temp.W=[.5 .5];
        
        X2.Values=zeros(X2.N,dim);
        [X2,Proposal_temp] =ISEM_astro(X2,Proposal_temp,Proposal0,df,Target_function,gama,j,data);
        
        
        
        ESS_temp=1/sum(X2.NormalizedWeight.^2)/X2.N;
        Sigma_temp=Proposal.Sigma;
        Sigma_temp(:,:,end+1:end+2)=Proposal_temp.Sigma;
        
        
        str1=num2str(j);
        str2='th iter--';
        str3='try to add a new component';
        strf = [str1 str2 str3 str0];
        if plot_id
            Plot_GM([X.Values;X2.Values],Proposal.M+1,[Proposal.W 0],[Proposal.Mu;Proposal_temp.Mu]',df/(df-2)*Sigma_temp(1:2,1:2,:));axis(axis_lmt);grid on; %text(text_Pos_X,text_Pos_Y,strf,'BackgroundColor',[.7 .9 .7]);
            aviobj = addframe(aviobj, getframe);
        end
        
        Proposal_try.M=Proposal.M+2;
        Proposal_try.Mu=[Proposal.Mu;Proposal_temp.Mu];
        Proposal_try.Sigma=Proposal.Sigma;
        Proposal_try.Sigma(:,:,Proposal_try.M-1:Proposal_try.M)=Proposal_temp.Sigma;
        Proposal_try.W=[(1-X2.N/(X.N+X2.N))*Proposal.W X2.N/(X.N+X2.N)*Proposal_temp.W];
        
        X_try.Values=[X.Values;X2.Values];
        X_try.N=X.N+X2.N;
        X_try.Resp=X.Resp;
        
        for i=1:Proposal.M
            X_try.Resp(X.N+1:X.N+X2.N,i)=exp(log_t_pdf(X2.Values,Proposal.Mu(i,:),Proposal.Sigma(:,:,i),df));
        end
        rp_temp(:,1)=exp(log_t_pdf(X.Values,Proposal_temp.Mu(1,:),Proposal_temp.Sigma(:,:,1),df));
        rp_temp(:,2)=exp(log_t_pdf(X.Values,Proposal_temp.Mu(2,:),Proposal_temp.Sigma(:,:,2),df));
        X_try.Resp=[X_try.Resp(1:X.N+X2.N,:) [rp_temp;X2.Resp]];
        
        index_survival=find(Proposal_try.W>w_thr);
        Proposal_try.M=length(index_survival);
        Proposal_try.Mu=Proposal_try.Mu(index_survival,:);
        Proposal_try.Sigma=Proposal_try.Sigma(:,:,index_survival);
        Proposal_try.W=Proposal_try.W(index_survival);
        Proposal_try.W=Proposal_try.W/sum(Proposal_try.W);
        
        X_try.Resp=X_try.Resp(:,index_survival);
        X_try.Proposal=X_try.Resp*Proposal_try.W';
        X_try.logProposal=log(X_try.Proposal);
        X_try.logAnnealTarget=[X.logAnnealTarget;X2.logAnnealTarget];
        X_try.logWeight=X_try.logAnnealTarget-X_try.logProposal;
        
        logWeight_scaled=X_try.logWeight-max(X_try.logWeight)+10; % Scaling
        weight_temp=exp(logWeight_scaled);
        X_try.NormalizedWeight=weight_temp/sum(weight_temp);
         
        Proposal=Proposal_try;
        X=X_try;
        ESS(j)=1/sum(X_try.NormalizedWeight.^2)/(X_try.N);
        
        
        str2='th iter--';
        str3=num2str(j2);
        str4='th new component added:ESS=';
        str5=' %i \r';
        strf = [str1 str2 str3 str4 str5];
        fprintf(strf, ESS(j));
        
        
        str2='th iter--';
        str3=num2str(j2);
        str4='th new component added';
        
        strf = [str1 str2 str3 str4 str0];
        if plot_id
            Plot_GM(X.Values(:,1:2),Proposal.M,Proposal.W,Proposal.Mu(:,1:2)',df/(df-2)*Proposal.Sigma(1:2,1:2,:));axis(axis_lmt);grid on; %text(text_Pos_X,text_Pos_Y,strf,'BackgroundColor',[.7 .9 .7]);
            aviobj = addframe(aviobj, getframe);
        end
        
        j2=j2+1;
        ESS_add=[ESS_add ESS(j)];
        if mod(j2,10)==0 && j2>=10
            tmptmp=sum(ESS_add(end-10+1:end-5));
            if (sum(ESS_add(end-5+1:end))-tmptmp)/tmptmp<.1
                break;
            end
        end
    end
    
    while 1
        ESS_em=ESS(j);
        X.N=N;
        X.Values=zeros(X.N,dim);
        X.Root=zeros(X.N,1);
        %Sample
        [X,Proposal]=ISEM_astro(X,Proposal,Proposal0,df,Target_function,gama,j,data);
        
        ESS(j)=1/sum(X.NormalizedWeight.^2)/(X.N);
        if ESS(j)<ESS_em
            break;
        end
    end
    
    %end
    str2='th iter--';
    str3=num2str(j2);
    str4='th new component added:ESS=';
    str5=' %i \r';
    strf = [str1 str2 str3 str4 str5];
    fprintf(strf, ESS(j));
    
    
    str2='th iter--';
    str3=num2str(j2);
    str4='th new component added';
    
    strf = [str1 str2 str3 str4 str0];
    if plot_id
        Plot_GM(X.Values(:,1:2),Proposal.M,Proposal.W,Proposal.Mu(:,1:2)',df/(df-2)*Proposal.Sigma(1:2,1:2,:));axis(axis_lmt);grid on; %text(text_Pos_X,text_Pos_Y,strf,'BackgroundColor',[.7 .9 .7]);
        aviobj = addframe(aviobj, getframe);
    end
    ESS_Arr=[ESS_Arr ESS(j)];
    
    fprintf('SAIS_final_rastrigin:j/Temp/ESS/M/N %i / %1.4f / %i /%i \r', j,gama(j),ESS(j),Proposal.M);
    
end

%-------------------importance sampling for computing marginal likelihood-------------%

Y.N=X.N;
Y=t_mix_sample(Proposal,Y,df);
Y.Resp=zeros(Y.N,Proposal.M); % responsibilily of each component with regard to each particle
for i=1:Proposal.M
    Y.Resp(:,i)=exp(log_t_pdf(Y.Values,Proposal.Mu(i,:),Proposal.Sigma(:,:,i),df));
end
Y.Proposal=Y.Resp*Proposal.W';
Y.logProposal=log(Y.Proposal);


[Y.logPrior,Y.logLike]=feval(Target_function,Y.Values,data);
Y.logTarget=Y.logPrior+Y.logLike;



Y.logWeight=Y.logTarget-Y.logProposal;
logWeight_scaled=Y.logWeight-max(Y.logWeight)+10; % Scaling
weight_temp=exp(logWeight_scaled);
Y.NormalizedWeight=weight_temp/sum(weight_temp);

ESS_final=1/sum(Y.NormalizedWeight.^2)/Y.N;
q=exp(Y.logWeight);
MarLik=mean(exp(Y.logWeight)) % marginal likelihood of this model

str1 = 'title(';
str2 = ['Particle Size N=' num2str(N)];
strf = [str1 'str2' ')'];
figure,plot(ESS);xlabel('Iteration');ylabel('ESS/N');
eval(strf);
