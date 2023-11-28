function Output = RosenbrockFunction(Samples)

Summation = zeros(size(Samples,1), 1);

for i = 1:size(Samples,2)-1

    Summation = Summation + (100 * (Samples(:, i + 1) - Samples(:, i) .^ 2) + (1 - Samples(:, i)) .^ 2);
end

Output = Summation;

end