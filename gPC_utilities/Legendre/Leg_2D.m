function [f] = Leg_2D(x1,x2,index);
% 
% Leg_2D     returns the value of the 2D Legendre polynomial of total order
% sum(index) evaluated at point(s) (x1,x2).
%
% Synopsis:  [f] = Leg_2D(x1,x2,index);
%
% Inputs:    x1,x2 = points abscissa
%            index = sequence of Legendre polynomial order in each
%            direction
% Output:    f = Legendre polynomial values at the point
%
f = Leg_1D_P15(x1,index(1)).*Leg_1D_P15(x2,index(2));