function vec = dec2base_(val, d, n)
vec(n) = 0;
for i = 0:n-1
  i_ = d^i;
  vec(end-i) = mod(floor(val / i_), d);
end
end
