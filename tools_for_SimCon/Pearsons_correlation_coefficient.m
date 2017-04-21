function [rho,pvalue] = Pearsons_correlation_coefficient(A,B)
%Computes the Pearson's correlation coefficient of 
% the vectorized upper triangles of input matrices A and B
%
% Author: Jason E. Hill, Ph.D.
% Date: 25 Jan 2017

N = size(A,1);
if N ~=size(A,2)
   error('Input matrices must be square and the same size!');
end
if N ~=size(B,1)
   error('Input matrices must be square and the same size!');
end
if N ~=size(B,2)
   error('Input matrices must be square and the same size!');
end

upperTriangle = triu(ones(N,N),1);

a = A(logical(upperTriangle));
b = B(logical(upperTriangle));

[r,p] = corrcoef(a,b);

rho = r(2);
pvalue = p(2);

end

