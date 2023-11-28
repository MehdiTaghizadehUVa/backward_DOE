function Y = IshigamiFunction(X)

if size(X, 2) ~= 3
    error('Ishigami model requires 3-dimensional inputs!');
end

a = 7;
b = 0.1;
X = X * pi;
Y = sin(X(:, 1)) + a * sin(X(:, 2) .^ 2) + (b * X(:, 3) .^ 4) .* sin(X(:, 1));

end