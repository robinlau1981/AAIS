
% Takes a parameter struct theta and a time vector t.  Returns the predicted radial velocities as a vector.

function radialVelocity = model_v4(theta, t)
% theta: [C K P e omega Mc0 K2 P2 e2 omega2 Mc02]
% compute mean anomaly: 

Mc = ((2*pi)/theta(3))*t + theta(6);
Mc2 = ((2*pi)/theta(3+5))*t + theta(6+5);

if (theta(4) < 0.001)&&(theta(4+5) < 0.001)
  radialVelocity = theta(1) + theta(2)*cos(Mc)+ theta(2+5)*cos(Mc2);
elseif (theta(4+5) < 0.001)
  % convert mean anomaly to eccentric anomaly:
  E = meananomaly2eccentricanomaly(Mc, theta(4));  
  % convert eccentric anomaly to true anomaly: 
  T = eccentricanomaly2trueanomaly(E,theta(4));  
  
  radialVelocity = theta(1) + theta(2)*(cos(theta(5)+T)+theta(4)*cos(theta(5)))+ theta(2+5)*cos(Mc2);
elseif (theta(4) < 0.001)
  % convert mean anomaly to eccentric anomaly:
  E2 = meananomaly2eccentricanomaly(Mc2, theta(4+5));
  
  % convert eccentric anomaly to true anomaly: 
  T2 = eccentricanomaly2trueanomaly(E2,theta(4+5));
  %T2=T2';
  radialVelocity = theta(1) + theta(2+5)*(cos(theta(5+5)+T2)+theta(4+5)*cos(theta(5+5)))+ theta(2)*cos(Mc);
else
  % convert mean anomaly to eccentric anomaly:
  E = meananomaly2eccentricanomaly(Mc, theta(4));
  
  % convert eccentric anomaly to true anomaly: 
  T = eccentricanomaly2trueanomaly(E,theta(4)); 
  %T=T';
  % convert mean anomaly to eccentric anomaly:
  E2 = meananomaly2eccentricanomaly(Mc2, theta(4+5));
  
  % convert eccentric anomaly to true anomaly: 
  T2 = eccentricanomaly2trueanomaly(E2,theta(4+5));
  %T2=T2';
  radialVelocity = theta(1)+ theta(2)*(cos(theta(5)+T)+theta(4)*cos(theta(5))) + theta(2+5)*(cos(theta(5+5)+T2)+theta(4+5)*cos(theta(5+5)));
end