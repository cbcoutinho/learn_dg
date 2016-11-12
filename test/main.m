N = 3;

figure; clf; hold on; %grid on;

x = linspace(-1, 1, 20);
for ii = 0:N
  plot(x, legendre(x, ii), 'o-')
%  plot(x, chebyshev(x, ii), 'o-')
end
