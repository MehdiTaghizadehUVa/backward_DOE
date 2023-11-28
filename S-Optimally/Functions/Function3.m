function Y = Function3(X)

dim = size(X, 2);
Y = zeros(size(X, 1), 1);

for i= 1:dim
    
  c = 1 / i^2; 
  Y = Y + c * (X(:, i) + 1) / 2 ;

end
pwr = -(dim +1);
Y = (1 + Y) .^ pwr;