% compute samle experience covariance
function Cov=Sample_experience_cov(x)
[r,l]=size(x);
Cov=zeros(r,r);
mx=mean(x,2);
for i=1:l
    Cov=Cov+(x(:,i)-mx)*(x(:,i)-mx)';
end
Cov=Cov./l;