function Y=Cmp_min_dis(Mu_candidate,Mu_cmps)
[r,c]=size(Mu_cmps);
for k=1:r
    dis(k)=norm(Mu_candidate-Mu_cmps(r,:));
end
[min_dis,min_ind]=min(dis);
Y=Mu_cmps(min_ind(1),:);