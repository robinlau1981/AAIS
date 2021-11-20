% takes a chain struct whose components are lists, and makes it a list whose components are theta structs.

function thetas = chains2thetas(chains)

n = length(chains.C);

for i = 1:n
  theta.C = chains.C(i);
  theta.K = chains.K(i);
  theta.P = chains.P(i);
  theta.e = chains.e(i);
  theta.omega = chains.omega(i);
  theta.M0 = chains.M0(i);
  theta.s = chains.s(i);
  theta.x = chains.x(i);
  theta.y = chains.y(i);
  theta.z = chains.z(i);
  theta.logP = chains.logP(i);
  theta.logK = chains.logK(i);
%  theta.logs = chains.logs(i);
  theta.l = chains.l(i);

  thetas(i) = theta;
end


