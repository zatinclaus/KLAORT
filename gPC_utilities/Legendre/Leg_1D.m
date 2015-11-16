function [value] = Leg_1D(degree, x)
% 
% Leg_1D    returns the value of the 1D Legendre polynomial of order P,
%             evaluated at point(s) x.
%
% Synopsis:  [value] = Leg_1D(degree, x);
%
% Inputs:    x = points abscissa
%            degree = Legendre polynomial order
% Output:    value = Legendre polynomial values at the point
%
% Remark:    explicit expressions up to P=10, then recurrence relation for
% higher order
%

switch(degree)
    case 0,
        value = ones(length(x),1);
    case 1,
        value = x;
    case 2,
        value = 0.5.*(3.*x.^2-1);
    case 3,
        value = 0.5.*(5.*x.^3-3.*x);
    case 4,
        value = 0.125.*(35.*x.^4 - 30.*x.^2 + 3);
    case 5,
        value = 0.125.*(63.*x.^5 - 70.*x.^3 + 15.*x);
    case 6,
        value = 0.0625.*(231.*x.^6 - 315.*x.^4 + 105.*x.^2 - 5);
    case 7,
        value = 0.0625.*(429.*x.^7 - 693.*x.^5 + 315.*x.^3 - 35.*x);
    case 8,
        value = (6435.*x.^8 - 12012.*x.^6 + 6930.*x.^4 -1260.*x.^2 + 35)./128;
    case 9,
        value = 12155/128.*x.^9 - 6435/32.*x.^7 + 9009/64.*x.^5 - 1155/32.*x.^3 + 315/128.*x;
    case 10,
        value = 46189/256.*x.^10 - 109395/256.*x.^8 + 45045/128.*x.^6 - 15015/128.*x.^4 +3465/256.*x.^2 - 63/256;
    otherwise
        d1 = degree-1; d2 = degree-2;
        value = ((2.*d1+1).*x.*Leg_1D(d1,x) - d1.*Leg_1D(d2,x))./degree;
end
