function n=func_PCE_NumPols(d,p);
%FUNC_PCE_NUMPOLS  Returns number of polynomials in full PCE
%
%       Syntax: n = func_PCE_NumPols(d,p); 
%
%       Input variables
%               d       dimension (integer)
%               p       order (integer)
%
%       Output variables
%               n       (integer)
%
%       Description
%               n=sum_{i=0}^{p}(i+d-1)!/i!/(d-1)!
%

%       Reference
%               C. Soize and R. Ghanem. Physical systems with random uncertainties: 
%               chaos representations with arbitrary probability measure. SIAM Journal 
%               on Scientific Computing, 26:395–410, 2004.
%
%		R. Ghanem and P. Spanos. Stochastic Finite Elements: A Spectral Approach. Springer, 1991.

%       Maarten Arnst, 01/21/2009
%       arnst@usc.edu


n=0;
for m1=0:p
   n=n+factorial(m1+d-1)/factorial(m1)/factorial(d-1);
end

