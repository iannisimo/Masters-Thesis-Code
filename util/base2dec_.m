function val = base2dec_(vec, d, n)
val = 0;
for i = 0:n-1
  i_ = d^i;
  val = val + (i_ * vec(end-i));
end
end
