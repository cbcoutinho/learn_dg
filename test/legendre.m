function L = legendre(x, N)

if N == 0
  L = ones(size(x));
elseif N == 1
  L = x;
else
  L = (2*N+1)/(N+1).*x.*legendre(x, N-1) - ... 
      N/(N+1).*legendre(x,N-2);
end
  
  
end