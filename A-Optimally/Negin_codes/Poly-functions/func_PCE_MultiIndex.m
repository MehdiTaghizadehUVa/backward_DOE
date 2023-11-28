function n=func_PCE_MultiIndex(j, d);
%FUNC_PCE_MULTIINDEX  Returns multi-index of multi-dimensional polynomial
%
%       Syntax: n = func_PCE_MultiIndex(j, d); 
%
%       Input variables
%               j       index (integer)
%               d       dimension (integer)
%
%       Output variables
%               n       coefficients (1*d integer vector)
%
%       Description
%               H_{j}(x)=h_{n(1)}(x(1))*...*h_{n(d)}(x(d))

%       Reference
%               C. Soize and R. Ghanem. Physical systems with random uncertainties: 
%               chaos representations with arbitrary probability measure. SIAM Journal 
%               on Scientific Computing, 26:395–410, 2004.
%
%               R. Ghanem and P. Spanos. Stochastic Finite Elements: A Spectral Approach. Springer, 1991.

%       Maarten Arnst, 01/21/2009
%       arnst@usc.edu

% Compute the polynomial order
index_p=0;
index_pp=1;
p=0;
while j>=index_pp,
    index_p=index_pp;
    index_pp=index_pp+factorial((p+1)+d-1)/factorial(p+1)/factorial(d-1);
    p=p+1;
end

% Compute the first index
index_pp=index_p+1;
newindex=0;
n=zeros(1,d);
for m1=0:p
    if j<index_pp
        n(1)=p-m1;
        if d~=1
            for m2=0:m1-1
                newindex=newindex+factorial(m2+d-2)/factorial(m2)/factorial(d-2);
                n(2:end)=func_PCE_MultiIndex(newindex+j-index_p, d-1);
            end
        end
        break;
    else
        index_p = index_pp;
        index_pp=index_pp+factorial(m1+1+d-2)/factorial(m1+1)/factorial(d-2);
    end
end
