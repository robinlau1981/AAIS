function pdf =log_Gaussian_pdf(X,Mu,Sigma) 
pdf=(normpdfln(X',Mu',[],Sigma))'; % Colomn vec
