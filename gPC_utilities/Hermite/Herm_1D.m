function [value] = Herm_1D(degree, x)
% 
% Herm_1D    returns the value of the 1D Hermite polynomial of order P,
%             evaluated at point(s) x.
%
% Synopsis:  [value] = Herm_1D(degree, x);
%
% Inputs:    x = points abscissa
%            P = Hermite polynomial order
% Output:    value = Hermite polynomial values at the point
%
% Remark:    explicit expressions up to P=10, then recurrence relation for
% higher order
%

switch(degree)
    case 0,
        % value = ones(length(x),1);
        value = 1+0.*x; % faster!
    case 1,
        value = x;
    case 2,
        value = x.^2-1;
    case 3,
        value = x.^3-3.*x;
    case 4,
        value = x.^4-6.*x.^2+3;
    case 5,
        value = x.^5-10*x.^3+15*x;
    case 6,
        value = x.^6-15*x.^4+45*x.^2-15;
    case 7,
        value = x.^7-21*x.^5+105*x.^3-105*x;
%     case 8,
%         value = x.^8-28*x.^6+210*x.^4-420*x.^2+105;
%     case 9,
%         value = x.^9-36*x.^7+378*x.^5-1260*x.^3+945*x;
%     case 10,
%         value = x.^10-45*x.^8+630*x.^6-3150*x.^4+4725*x.^2-945;        
    otherwise
        value = x.*Herm_1D(degree-1,x) - (degree-1).*Herm_1D(degree-2,x);
end
