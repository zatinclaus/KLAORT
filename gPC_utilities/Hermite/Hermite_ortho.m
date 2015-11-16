function [value] = Hermite_ortho(N, p, mi, x)
%
% Hermite    returns the value of the N-dimensional ORTHONORMAL Hermite polynomial of order p,
%             evaluated at point(s) x.
%
% Synopsis:  [value] = Hermite_ortho(N, p, mi, x);
%
% Inputs:    N = stochastic dimension
%            mi = basis multi-indices [P,N]
%            x = points abscissa [Nquad,N]
%            p = Hermite polynomial order (from 0 to P)
% Output:    value = Hermite orthonormal polynomial values at the points
%
% Remark:
%

if N == 1,
    value = Herm_1D(p,x);
    value = value./sqrt((factorial(mi(p+1))));
else
    x_size = size(x);
    temp = zeros(x_size(1),N);
    for j=1:N,
        temp(:,j) = Herm_1D(mi(p+1,j),x(:,j));
    end
    value = prod(temp,2);
    % value = value./sqrt(prod(factorial(mi(p+1,:))));
    value = value./sqrt(prod(my_factorial(mi(p+1,:)))); % slightly faster version
end

end
