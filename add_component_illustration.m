clear;
close all;
N=1000;

a=normrnd(20,4,1,N);
[fa,xa] = ksdensity(a);

b=normrnd(0,4,1,N);
[fb,xb] = ksdensity(b);

for i=1:N
    if rand<.5
        c(i)=normrnd(0,4);
    else
        c(i)=normrnd(20,4);
    end
end
[fc,xc] = ksdensity(c);
figure,
plot(xa,fa,'k',xb,fb,'r',xc,fc,'b');legend('target distribution','original proposal','proposal after adaptation')%axis([-50 50 0 .1]);
hold on;plot([0 0],[0 .11],[20 20],[0,.09]);hold off;

gtext('the new added component');gtext('the location of the sample with biggest importance weight');gtext('the location of the sample with smallest importance weight');
