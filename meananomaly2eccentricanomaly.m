% I think this is the best version yet of this function.  Any e works, any
% M or vector of Ms.

function E = meananomaly2eccentricanomaly(M, e)

E = zeros(length(M),1);

decimalplaceaccuracy = 8;
% testing showed that you hardly lost any speed by going from 3 decimal place accuracy to 8 decimal place accuracy.  So I've made this a constant in the function rather than an argument.

M_is_between_pi_and_2pi = 0;

if e == 0
    E = mod(M,2*pi);
else
    for i = 1:length(M)
        Mtemp = M(i);
        Mtemp = mod(Mtemp,2*pi);
        if Mtemp > pi
            Mtemp = 2*pi-Mtemp;
            M_is_between_pi_and_2pi = 1;
        end
        if e>0.9 & Mtemp~=0    % this is to mend a situation when e near 1 and Mtemp near 0 results in an 'offshooting'.  This is a better seed.
            x0 = exp(0.67321563757256+0.38477427304553*log(Mtemp));
        else
            x0 = Mtemp;
        end
        %x0 = Mtemp;
        count = 1;
        while abs(Mtemp - x0 + e*sin(x0)) > (10^-decimalplaceaccuracy) % | isnan(x0)
            count = count + 1;
            if count > 5
                x0 = rand*pi;
                count = 1;
                %disp(['had to reboot for e = ',num2str(e),', M = ',num2str(M(i)),', Mtemp = ',num2str(Mtemp),'.'])
            end

            A = 0.5*e*sin(x0);
            B = 1 - e*cos(x0) - x0*e*sin(x0);
            C = -Mtemp + x0 - e*sin(x0) -x0*(1 -e*cos(x0)) + 0.5*(x0^2)*(e*sin(x0));
            x0 = -B/(2*A) + sqrt(B^2 - 4*A*C)/(2*A);

        end % while
        if M_is_between_pi_and_2pi
            E(i) = 2*pi-x0;
            M_is_between_pi_and_2pi = 0;
        else
            E(i) = x0;
        end
    end % for loop
end %  if e==0
