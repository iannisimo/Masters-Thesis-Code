function [ctrl, ctrlVal, targ, V] = singleClubHouseholder(term, psi_j, d)
  term = char(term);
  ctrl = -1;
  ctrlVal = 0;
  targ = find(term == '-', 1);
  ctrl_ = term > '0';
  if any(ctrl_)
    ctrl = find(ctrl_, 1, 'last');
    ctrlVal = term(ctrl);
  end
  ctrlterm = term(term ~= '-');
  phi = zeros(size(psi_j));
  for k = 0:d-1
    t_ = [ctrlterm, char('0' + k), repmat('0', 1, length(term) - length(ctrlterm) - 1)];
    l = zeros(size(psi_j));
    l((t_ - '0') * d.^(length(t_)-1:-1:0)' + 1, 1) = 1;
    r = zeros(size(psi_j));
    r(k+1, 1) = 1;
    phi = phi + l' * psi_j * r;
  end
  V = makeHouseholder(phi(1:d));
end