function Output = LD_HOFunction(Samples)

Summation = zeros(size(Samples,1), 1);

for i = 1:5

Summation = Summation + (Samples(:, 1) .^ (2 * i)) .* (Samples(:, 2) .^ (2 * i));
end

Output = Summation;
end