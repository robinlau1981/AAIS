function mvec=rv_1p(theta,t)
nd=length(data.V);
aw1=.5*(theta(4)-theta(6));
fvec=mod(t/theta(1)+aw1,1);
global a2;
a2=theta(2);

q_start=0;
q_end=1;
X0=0;
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[q,X]= ode45(@transfun,[q_start,q_end],X0,options);

Xvec=zeros(nd,1);
for k=1:nd
    qc=fvec(k); % candidate value of theta
    ind1=find(q==qc);
    if isempty(ind1)~=1
        Xvec(k)=X(ind1);
    else
        qs=sort(q);
        ind2=find(qs>qc);
        Xvec(k)=X(min(ind2)-1)+(qc-q(min(ind2)-1))*(X(min(ind2))-X(min(ind2)-1))/(q(min(ind2))-q(min(ind2)-1));
    end
end

mvec=theta(3)+theta(5)*cos((Xvec+aw1));


function dtheta_dq=transfun(q,X)
% See Gregory3plan47UMaParallel21Feb10M7.nb for definition
global a2;
dtheta_dq=2*pi/(1-a2^2)^1.5*(1+a2*cos(X))^2;