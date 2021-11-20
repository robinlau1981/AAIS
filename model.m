
% Takes a parameter struct theta and a time vector t.  Returns the predicted radial velocities as a vector.

function radialVelocity = model(theta, t)

% compute mean anomaly: 

M = ((2*pi)/theta.P)*t + theta.M0;

if theta.e < 0.001
  radialVelocity = theta.C + theta.K*cos(M);
else
  % convert mean anomaly to eccentric anomaly:
  E = meananomaly2eccentricanomaly(M, theta.e);
  % convert eccentric anomaly to true anomaly: 
  T = eccentricanomaly2trueanomaly(E,theta.e);
  radialVelocity = theta.C + theta.K*(cos(theta.omega+T)+theta.e*cos(theta.omega));
end



















