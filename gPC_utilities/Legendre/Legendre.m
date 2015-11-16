function [value] = Legendre(N, p, mi, x)
%
% Hermite    returns the value of the N-dimensional Legendre polynomial of order p,
%             evaluated at point(s) x.
%
% Synopsis:  [value] = Legendre(N, p, mi, x);
%
% Inputs:    N = stochastic dimension
%            mi = basis multi-indices [P,N]
%            x = points abscissa [Nquad,N]
%            p = Legendre polynomial order (from 0 to P)
% Output:    value = Legendre polynomial values at the points
%
% Remark:
%

if N == 1,
    value = Leg_1D(p,x);
else
    x_size = size(x);
    temp = zeros(x_size(1),N);
    for j=1:N,
        temp(:,j) = Leg_1D(mi(p+1,j),x(:,j));
    end
    value = prod(temp,2);
end


end
