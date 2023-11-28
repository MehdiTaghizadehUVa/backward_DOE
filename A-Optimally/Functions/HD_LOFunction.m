function Output = HD_LOFunction(Samples)

Summation = zeros(size(Samples,1), 1);

for i = 1:size(Samples,2)-1

Summation = Summation + Samples(:, i) .* Samples(:, i+1);
end

Output = Summation;
end