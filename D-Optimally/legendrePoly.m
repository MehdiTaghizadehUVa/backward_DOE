function A = legendrePoly(order, X)
number_Sample = size(X, 1);
Pol_dim = size(X, 2);
Pol_Type = "Legendre";
% number of coeffient
n_coef = factorial(Pol_dim + order) / (factorial(Pol_dim) * factorial(order));
A = ones(number_Sample, n_coef);        
%construct the measurement matrix
for i=1:n_coef
        multi_index = func_PCE_MultiIndex(i-1, Pol_dim); 
            for d=1:Pol_dim
            c=func_PCE_1DPoly(multi_index(d), Pol_Type);
                for j=1:number_Sample
                A(j,i) = polyval(c,X(j,d))*A(j,i);
                end      
            end
        end

end