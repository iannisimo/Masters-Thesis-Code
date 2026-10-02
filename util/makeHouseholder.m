function [W] = makeHouseholder(phi)
zero = zeros(size(phi));
zero(1,1) = 1;
dot = zero'*phi;
if abs(dot) >= (1 - eps) * norm(phi)
  % phi is zero or already a multiple of e1: nothing to reflect
  W = eye(length(phi));
  return
end
if abs(dot) < eps * norm(phi)
  sgn = 1;
else
  sgn = dot / abs(dot);
end
% phi + ||phi|| sgn e1 avoids the cancellation in eta(1) when phi is close
% to e1; the reflector then maps phi to -sgn ||phi|| e1, so it is negated
% to keep W phi = sgn ||phi|| e1
eta = phi + sqrt(phi'*phi) * sgn * zero;
W = -(eye(length(phi)) - (2 / (eta'*eta)) * (eta * eta'));
end
