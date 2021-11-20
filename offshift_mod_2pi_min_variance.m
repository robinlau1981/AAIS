
% returns the amount off offshift required to minimize the variance of a list that is mod 2*pi.

function offshift = offshift_mod_2pi_min_variance(phis)

offshift_candidates = [0:0.01:2*pi];
variances = zeros(length(offshift_candidates),1);
for i = 1:length(offshift_candidates)
  variances(i) = var(mod(phis+offshift_candidates(i),2*pi));
end

[smallest, index] = min(variances);

offshift = offshift_candidates(index);



