function Y = SobolFunction(X)
a = [1 2 5 10 20 50 100 500];
X = (X + 1) / 2;
dim = size(X, 2);
Y = ones(size(X, 1), 1);
for i= 1:dim
  Y = Y .* ((abs(4 * X(:, i) - 2) + a(i)) / (1 + a(i))) ;
end