N=1e3;
e=.2;
M=unifrnd(0,2*pi,1,N);
tic;
for k=1:N
E(k)=fzero(@(E) trans_eqn(E,e,M(k)),0);
end
time1=toc
tic;
E2 = meananomaly2eccentricanomaly(M, e);
time2=toc