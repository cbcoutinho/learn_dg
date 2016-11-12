function T = chebyshev(x, N)
  
  
if N == 0
  T = ones(size(x));
else if N == 1
  T = x;
else
  T = 2.*x.*chebyshev(x, N-1) - chebyshev(x, N-2);
end

end