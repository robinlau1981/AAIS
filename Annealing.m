function Ti=Annealing(T0,Tn,N,iter, schedule);
% get the temperature of defferent annealing stratgies.
%Refer to Brian T. Luke, Simulated Annealing Cooling Schedules,
%http://fconyx.ncifcrf.gov/~lukeb/simanf1.html.
switch schedule
    case 1
        Ti=T0-iter*(T0-Tn)/N;
    case 2
        Ti=T0*(Tn/T0).^(iter/N);
    case 3
        A=(T0-Tn)*(N+1)/N;
        B=T0-A;
        Ti=A./(iter./100+1)+B;
    case 4
        Ti=0.5*(T0-Tn)*(1+cos(iter*pi/N))+Tn;
    case 5
        Ti=0.5*(T0-Tn)*(1-tanh(10*iter/N-5))+Tn;
    case 6
        Ti=(T0-Tn)/cosh(10*iter/N)+Tn;
    case 7
        A=(1/(N*N))*log(T0/Tn);
        Ti=T0*exp(-A.*iter.^2);
%     case 8
%         A=2.5*(1/N)*log(T0/Tn);
%         Ti=T0*exp(-A*iter);
end
