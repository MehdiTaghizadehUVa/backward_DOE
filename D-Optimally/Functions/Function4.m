function Y = Function4(X)

dim = size(X, 2);
Y = ones(size(X, 1), 1);

c = [-3 , 2];
w = [0.5, 0,5];
for i= 1:dim
  Y = Y .* ((((X(:, i) + 1) / 2 - w(i)) .^ 2 + c(i) ^ -2) .^ -1);
end
end