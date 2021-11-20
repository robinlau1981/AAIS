function Xminus=unif_resample(X,q)
[N,dim]=size(X);
rand_num = rand; % uniform random number between 0 and 1
qtempsum = 0;
for ii = 1 : N
    qtempsum = qtempsum + q(ii);
    if qtempsum >= rand_num
        Xminus=X(ii,:);
        break;
    end
end