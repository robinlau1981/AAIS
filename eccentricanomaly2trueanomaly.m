% converts from eccentricanomaly to true anomaly.

function trueanomaly = eccentricanomaly2trueanomaly(E,e)

% Used to diagnose why I'm getting imaginary errors:
% if ~isreal(2*atan(tan(E/2)*sqrt((1+e)/(1-e))))
%   E
%   e
% end


trueanomaly = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));




