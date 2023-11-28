function Y = Function2(X)

dim = size(X, 2);
Y = zeros(size(X, 1), 1);

for i= 1:dim
  
  Y = Y + i * abs((X(:, i) + 1) / 2 - 0.5);

end

Y = exp(-Y);