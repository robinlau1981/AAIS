% compute samle experience covariance
function [mx,Cov]=Sample_experience_cov_v2(x)
x=x';
[r,l]=size(x);
Cov=zeros(r,r);
mx=mean(x,2);
for i=1:l
    Cov=Cov+(x(:,i)-mx)*(x(:,i)-mx)';
end
Cov=Cov./l;
mx=mx';