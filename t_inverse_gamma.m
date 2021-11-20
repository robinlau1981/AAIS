clear;
ndraws = 10000;
Lcut = 0.05;
Hcut = 0.95;
%
seednumber = 0; %4720;    %472534;   % if 0, random state at each clock time
           % good one 420 for [29 45], [29 54]
if seednumber
   randn('state',seednumber);
   rand('state',seednumber);
else
   randn('state',fix(100*sum(clock)));
   rand('state',fix(100*sum(clock)));
end
% %-----------------------------------------------------------------------------------
% %---------------------------- Gamma distribution ----------------------------------%
% %--- p(x) = ( b^a/Gamma(a) ) x^(a-1) exp(-bx) for a>0 and b>0.
% %---    where a is shape and b is inverse scale (rate) parameter.
% %--- E(x) = a/b;  var(x) = a/b^2;
% %--- Noninformative distribution: a,b -> 0.
% %--- The density function is finite if a >= 1.
% %-----------------------------------------------------------------------------------
% xm = 1.0;  %Mean
% xs = 0.5;  %Std
% 
% a = xm^2.0/xs^2.0;
% b = xm/xs^2.0;
% if (a<0.0) || (b<0.0)
%    error('plot_priordensity.m: the values of a and b implied by xm and xs are incorrect!');
% end
% 
% %a
% %b
% %[xm xs]
% %[a/b sqrt(a/b^2)]
% 
% x = 0.0:0.1:100.0;
% y = exp( a*log(b) - gammaln(a) + (a-1)*log(x) - b*x );
% %
% ydraws = gamrnd(a,1/b,[ndraws 1]);  %Note that in Matlab, it is 1/b, NOT b.
% ydraws_sort = sort(ydraws);
% yL = ydraws_sort(Lcut*ndraws);
% yH = ydraws_sort(Hcut*ndraws);
% 
% figure
% plot(x,y);
% str1 = 'title(';
% str2 = ['Gamma pdf of x with .90 interval (' num2str(yL) ',  ' num2str(yH) ')'];
% strf = [str1 'str2' ')'];
% eval(strf);
% %
% str1 = 'xlabel(';
% str2 = ['xm=' num2str(xm) ', xs=' num2str(xs) ' with a=' num2str(a) ' & b=' num2str(b)]; % ' and with .90 interval (' num2str(yL) '  ' num2str(yH) ')'];
% strf = [str1 'str2' ')'];
% eval(strf);



%-----------------------------------------------------------------------------------
%------------------------ Inverse-Gamma distribution ------------------------------%
%--- p(x) = ( b^a/Gamma(a) ) x^(-a-1) exp(-b/x) for a>0 and b>0.
%---    where a is shape and b is scale parameter.
%--- E(x) = b/(a-1) for a>1;  var(x) = b^2/( (a-1)^2*(a-2) ) for a>2;
%--- Noninformative distribution: a,b -> 0.
%--- How to draw: (1) draw z from Gamma(a,b); (2) let x=1/z.
%-----------------------------------------------------------------------------------
% xm = 30;  %Mean
% xs = 1;  %Std
% 
% a = xm^2.0/xs^2.0 + 2.0;
% b = (a-1.0)*xm;
% if (a<2.0) || (b<0.0)
%    error('plot_priordensity.m: the values of a and b implied by xm and xs are incorrect!');
% end
a=2; b=2;



%a=0.65;
%b=0.15;
%xm = b/(a-1);
%xs = sqrt(b^2/((a-1)^2*(a-2)));


x = 0.0:1e-3:10.0;
y = exp( a*log(b) - gammaln(a) - (a+1)*log(x) - b./x );
%
zdraws = gamrnd(a,1/b,[ndraws 1]);  %Note that in Matlab, it is 1/b, NOT b.
ydraws = 1.0./zdraws;
ydraws_sort = sort(ydraws);
yL = ydraws_sort(Lcut*ndraws);
yH = ydraws_sort(Hcut*ndraws);

figure
plot(x,y);
str1 = 'title(';
str2 = ['IG pdf of x with .90 interval (' num2str(yL) ',  ' num2str(yH) ')'];
strf = [str1 'str2' ')'];
eval(strf);
%
str1 = 'xlabel(';
str2 = ['xm=' num2str(xm) ', xs=' num2str(xs) ' with a=' num2str(a) ' & b=' num2str(b)]; % ' and with .90 interval (' num2str(yL) '  ' num2str(yH) ')'];
strf = [str1 'str2' ')'];
eval(strf);