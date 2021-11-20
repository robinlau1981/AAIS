
% March 20 2007: I'm adding a factor of K to this prior, which I think
% might correct for the fact that when I actually draw from the prior, I'm
% drawing *logK*, not K.  So the prior is of the transformed variables,
% includign logK.  I'm going to see if that makes the prior integrate to 1.


function l = log_prior_transformed_v5(thetas)

[r,c]=size(thetas);

% Y(1,j)=X(1,j);      % C
% Y(2,j)=exp(X(2,j)); % K
% Y(3,j)=exp(X(3,j)); % P
% Y(4,j)=sqrt(X(4,j)^2+X(5,j)^2); % e
% Y(5,j)=mod(atan(X(5,j)/X(4,j))+pi*(X(4,j)<0),2*pi); % omega
% Y(6,j)=mod(X(6,j)- Y(5,j),2*pi); % M0
% Y(12,j)=exp(X(7,j)); % s

P_min = 1;
P_max =2e4;% 365*1000;%1e3;%
K_0 = 1;
K_max = 2128; % Eric has a reason for this.
C_min = -K_max;
C_max = K_max;
s_0 = 1;
s_max = K_max;
log_kappa = - log(C_max-C_min) - 2*log(log(P_max/P_min)) - 2*log(log(1+(K_max/K_0))) - 2*log(K_0) - log(log(1+(s_max/s_0))) - log(s_0) - 2*2*log(2*pi);


l=zeros(r,1);

% for i = 1:r
% 
%     theta = thetas(i,:);
%     % conditions for falling in the right range; the indicator function for the prior:
%     if (C_min <= theta(1)) && (theta(1) <= C_max) && (0 < exp(theta(2))) && (theta(2) <= log(K_max)) && (P_min <= exp(theta(3))) &&...
%             (exp(theta(3)) <= P_max) && (theta(4)^2 + theta(5)^2 < 1) && (0 <= theta(6)) && (theta(6) < 2*pi)&& (0 < exp(theta(7))) && ...
%             (theta(7) <= log(K_max)) && (P_min <= exp(theta(8))) && (exp(theta(8)) <= P_max) && (theta(9)^2 + theta(10)^2 < 1) && (0 <= theta(11))...
%             && (theta(11) < 2*pi) && (0 < exp(theta(12))) && (exp(theta(12)) < s_max) && (theta(3)<theta(8))
%         l(i) = log_kappa + theta(2) - log(1+(exp(theta(2))/K_0))-0.5*log(theta(4)^2 + theta(5)^2)  + theta(7) - log(1+(exp(theta(7))/K_0))-0.5*log(theta(9)^2 + theta(10)^2)+ log(exp(theta(12))) - log(1+(exp(theta(12))/s_0)) ;
% 
%     else
% 
%         l(i) = -10^100;
%     end
% 
% end % for loop

for i = 1:r
    
% P_min = 1;
% P_max = 365*1000;%1e3;%
% K_0 = 1;
% K_max = 2128; % Eric has a reason for this.
% C_min = -K_max;
% C_max = K_max;
% s_0 = 1;
% s_max = K_max;
% log_kappa = - log(C_max-C_min) - 2*log(log(P_max/P_min)) - 2*log(log(1+(K_max/K_0))) - 2*log(K_0) - log(log(1+(s_max/s_0))) - log(s_0) - 2*2*log(2*pi);
%     
    
    
    % conditions for falling in the right range; the indicator function for the prior:
    if (C_min <= thetas(i,1)) && (thetas(i,1) <= C_max) && (0 < exp(thetas(i,2))) && (thetas(i,2) <= log(K_max)) && (P_min <= exp(thetas(i,3))) &&...
            (exp(thetas(i,3)) <= 2e4) && (thetas(i,4)^2 + thetas(i,5)^2 < 1) && (0 <= thetas(i,6)) && (thetas(i,6) < 2*pi)&& (0 < exp(thetas(i,7))) && ...
            (thetas(i,7) <= log(K_max)) && (P_min <= exp(thetas(i,8))) && (exp(thetas(i,8)) <= 2e4) && (thetas(i,9)^2 + thetas(i,10)^2 < 1) && (0 <= thetas(i,11))...
            && (thetas(i,11) < 2*pi) && (0 < exp(thetas(i,12))) && (exp(thetas(i,12)) < s_max) && (0 < exp(thetas(i,2+5*2))) && (thetas(i,2+5*2) <= log(K_max)) && (P_min <= exp(thetas(i,3+5*2))) &&...
            (exp(thetas(i,3+5*2)) <= 2e4) && (thetas(i,4+5*2)^2 + thetas(i,5+5*2)^2 < 1) && (0 <= thetas(i,6+5*2)) && (thetas(i,6+5*2) < 2*pi)
        l(i) = log_kappa + thetas(i,2) - log(1+(exp(thetas(i,2))/K_0))-0.5*log(thetas(i,4)^2 + thetas(i,5)^2)  + thetas(i,7) - log(1+(exp(thetas(i,7))/K_0))-0.5*log(thetas(i,9)^2 + thetas(i,10)^2)+ thetas(i,7+5) - log(1+(exp(thetas(i,7+5))/K_0))-0.5*log(thetas(i,9+5)^2 + thetas(i,10+5)^2)+ log(exp(thetas(i,17))) - log(1+(exp(thetas(i,17))/s_0)) ;

    else

        l(i) = -Inf;
    end

end % for loop


