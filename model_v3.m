% Takes a parameter struct theta and a time vector t.  Returns the predicted radial velocities as a vector.

function radialVelocity = model_v3(theta, t)

% compute mean anomaly: 

M = ((2*pi)/theta.P)*t + theta.M0;
M2 = ((2*pi)/theta.P2)*t + theta.M02;

if (theta.e < 0.001)&(theta.e2 < 0.001)
  radialVelocity = theta.C + theta.K*cos(M)+ theta.K2*cos(M2);
elseif theta.e2 < 0.001
  % convert mean anomaly to eccentric anomaly:
  E = meananomaly2eccentricanomaly(M, theta.e);
  % convert eccentric anomaly to true anomaly: 
  T = eccentricanomaly2trueanomaly(E,theta.e);
  radialVelocity = theta.C + theta.K*(cos(theta.omega+T)+theta.e*cos(theta.omega))+ theta.K2*cos(M2);
elseif theta.e< 0.001
  % convert mean anomaly to eccentric anomaly:
  E2 = meananomaly2eccentricanomaly(M2, theta.e2);
  % convert eccentric anomaly to true anomaly: 
  T2 = eccentricanomaly2trueanomaly(E2,theta.e2);
  radialVelocity = theta.C + theta.K2*(cos(theta.omega2+T2)+theta.e2*cos(theta.omega2))+ theta.K*cos(M);    
else
  E = meananomaly2eccentricanomaly(M, theta.e);
  % convert eccentric anomaly to true anomaly: 
  T = eccentricanomaly2trueanomaly(E,theta.e);    
  % convert mean anomaly to eccentric anomaly:
  E2 = meananomaly2eccentricanomaly(M2, theta.e2);
  % convert eccentric anomaly to true anomaly: 
  T2 = eccentricanomaly2trueanomaly(E2,theta.e2);
  radialVelocity = theta.C + theta.K*(cos(theta.omega+T)+theta.e*cos(theta.omega))+ theta.K2*(cos(theta.omega2+T2)+theta.e2*cos(theta.omega2)); 
end