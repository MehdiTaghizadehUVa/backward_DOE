function Output = PDEFunction(Samples)
number_Sample = size(Samples, 1);
Output = zeros(number_Sample, 1);

for i=1:number_Sample
    R = Samples(i,:);
    Output(i) = ODE(R);
end
end