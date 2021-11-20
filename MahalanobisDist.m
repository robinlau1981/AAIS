function D=MahalanobisDist(X,mu,Sigma)
[N,dim]=size(X);
D=zeros(1,N);
for i=1:N
    %D(i)=(X(i,:)-mu)*inv(Sigma)*(X(i,:)-mu)';
    D(i)=sqdist(X(i,:)',mu',inv_posdef(Sigma));
end