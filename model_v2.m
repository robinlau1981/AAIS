
% Takes a parameter struct theta and a time vector t.  Returns the predicted radial velocities as a vector.

function radialVelocity = model_v2(theta, t)
% theta: [C K P e omega M0]
% compute mean anomaly: 

M = ((2*pi)/theta(3))*t + theta(6);

if theta(4) < 0.001
  radialVelocity = theta(1) + theta(2)*cos(M);
else
  % convert mean anomaly to eccentric anomaly:
  E = meananomaly2eccentricanomaly(M, theta(4));  
  % convert eccentric anomaly to true anomaly: 
  T = eccentricanomaly2trueanomaly(E,theta(4)); 
  radialVelocity = theta(1) + theta(2)*(cos(theta(5)+T)+theta(4)*cos(theta(5)));
end



















