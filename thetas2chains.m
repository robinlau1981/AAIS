
% converts a list of theta structs into a struct of chains.

function chains = thetas2chains(thetas)

n = length(thetas);

for i = 1:n
  chains.C(i) = thetas(i).C;
  chains.K(i) = thetas(i).K;
  chains.P(i) = thetas(i).P;
  chains.e(i) = thetas(i).e;
  chains.omega(i) = thetas(i).omega;
  chains.M0(i) = thetas(i).M0;
  chains.s(i) = thetas(i).s;
  chains.logP(i) = thetas(i).logP;
  chains.logK(i) = thetas(i).logK;
  chains.logs(i) = thetas(i).logs;
  chains.x(i) = thetas(i).x;
  chains.y(i) = thetas(i).y;
  chains.z(i) = thetas(i).z;
  chains.l(i) = thetas(i).l;
end
