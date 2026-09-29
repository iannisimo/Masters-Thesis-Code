function seq = makeClubSequence(d, n)
  if n == 1, seq = "-"; return; end
  seq = [];
  seq_ = makeClubSequence(d, n-1);
  for q=0:d-1
    seq = [seq, seq_.insertBefore(1, char('0' + q))];
  end
  seq = [seq, string(repmat('-', 1, n))];
end