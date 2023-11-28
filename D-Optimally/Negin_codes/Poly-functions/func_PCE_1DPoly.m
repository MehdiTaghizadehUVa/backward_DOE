function c=func_PCE_1DPoly(p, varargin);
%       Returns coefficients of 1D normalized polynomial from Wiener-Askey scheme
%
%	Syntax: c = func_PCE_1DPoly(p, 'Hermite');
%		c = func_PCE_1DPoly(p, 'Legendre');
%		c = func_PCE_1DPoly(p, 'Laguerre', a);
%		c = func_PCE_1DPoly(p, 'Chebyshev');
%		c = func_PCE_1DPoly(p, 'Jacobi', a, b);  	
%
%	Input variables
%		p	degree (integer)
%               a       Laguerre parameter (real scalar > -1)
%               a, b    Jacobi parameters (real scalars > -1)
%
%	Output variables
%		c	coefficients (1*(p+1) real vector)
%
%	Description
%		Hermite orthogonal w.r.t. 	(2\pi)^{-1/2}\exp(-x^{2}/2)	on \real
%		Legendre orthogonal w.r.t. 	1				on ]-1,1[
%		Laguerre orthogonal w.r.t. 	x^{a}\exp(-x)			on ]0,+\infty[
%		Chebyshev orthogonal w.r.t.	(1-x^{2})^{-1/2}		on ]-1,1[
%		Jacobi orthogonal w.r.t.	(1-x)^{a}(1+x)^{b}		on ]-1,1[	  
%
%	See also
%		polyval

%	Reference
%               C. Soize and R. Ghanem. Physical systems with random uncertainties: 
%               chaos representations with arbitrary probability measure. SIAM Journal 
%               on Scientific Computing, 26:395–410, 2004.
%
%               R. Ghanem and P. Spanos. Stochastic Finite Elements: A Spectral Approach. Springer, 1991.

%	Maarten Arnst, 01/21/2009
%	arnst@usc.edu
if strcmp(varargin{1},'Hermite')

    switch p
    case 0,
        c=[1];
    case 1,
        c=[1 0];
    otherwise,
    % Higher order polynomials using recursion formula
        He_n=[1 0];
        He_nminus1=[1];
        for m1=2:p
            c=conv([1 0],He_n)-(m1-1)*[zeros(1,2) He_nminus1];
            if m1~=p
                He_nminus1=He_n;
                He_n=c;
            end
        end
    end
    c=(1/sqrt(factorial(p)))*c; % normalization
elseif strcmp(varargin{1},'Legendre')
    switch p
    case 0,
        c=[1];
    case 1,
        c=[1 0];
    otherwise
    % Higher order polynomials using recursion formula
        He_n=[1 0];
        He_nminus1=[1];
        for m1=2:p
            c=((2*m1-1)/m1)*conv([1 0],He_n)-((m1-1)/m1)*[zeros(1,2) He_nminus1];
            if m1~=p
                He_nminus1=He_n;
                He_n=c;
            end
        end
    end

    c=sqrt((2*p+1)/2)*c; % normalization
elseif strcmp(varargin{1},'Laguerre')
    a=varargin{2};
    switch p
        case 0,
            c=[1];
        case 1,
            c=[-1 a+1];
    otherwise
    % Higher order polynomials using recursion formula
        He_n=[-1 a+1];
        He_nminus1=[1];
        for m1=2:p
            c=(-1/m1)*conv([1 0],He_n)+((2*m1-1+a)/m1)*[0 He_n]-((m1-1+a)/m1)*[zeros(1,2) He_nminus1];
            if m1~=p
                He_nminus1=He_n;
                He_n=c;
            end
        end
    end
    c=sqrt(factorial(p)/gamma(a+p+1))*c; % normalization
elseif strcmp(varargin{1},'Chebyshev')
    switch p
    case 0,
        c=[1];
    case 1,
        c=[1 0];
    otherwise
    % Higher order polynomials using recursion formula
        He_n=[1 0];
        He_nminus1=[1];
        for m1=2:p
            c=2*conv([1 0],He_n)-[zeros(1,2) He_nminus1];
            if m1~=p
                He_nminus1=He_n;
                He_n=c;
            end
        end
    end
    if p==0
        c=sqrt(1/pi)*c; % normalization
    else
        c=sqrt(2/pi)*c;
    end
elseif strcmp(varargin{1},'Jacobi')
    a=varargin{2};b=varargin{3};
    switch p
    case 0,
        c=[1];
    case 1,
        c=[(a+b)/2+1 (a-b)/2];
    otherwise
        % Higher order polynomials using recursion formula
        He_n=[(a+b)/2+1 (a-b)/2];
        He_nminus1=[1];
        for m1=2:p
            alpha=(a+b+2*m1-1)*((a+b+2*m1-2)*(a+b+2*m1))/(2*m1*(a+b+m1)*(a+b+2*m1-2));
            beta=(a+b+2*m1-1)*(a^2-b^2)/(2*m1*(a+b+m1)*(a+b+2*m1-2));
            gammaa=2*(a+m1-1)*(b+m1-1)*(a+b+2*m1)/(2*m1*(a+b+m1)*(a+b+2*m1-2));
            c=alpha*conv([1 0],He_n)+beta*[0 He_n]-gammaa*[zeros(1,2) He_nminus1];
            if m1~=p
                He_nminus1=He_n;
                He_n=c;
            end
        end
    end
    c=sqrt(factorial(p)*(a+b+1+2*p)*gamma(a+b+p+1)/2^(a+b+1)/gamma(a+p+1)/gamma(b+p+1))*c; % normalization
end


