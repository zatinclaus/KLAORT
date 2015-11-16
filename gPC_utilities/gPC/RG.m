function [Xi] = RG(N,n,dist_type);
% this function generates n dimensions of N random numbers according to matlab
% random number generators for different distributions dist_type
% dist_type = 1:  % Normal Gaussian distribution
% dist_type = 2:  % Uniform distribution

switch dist_type,
case 1
  Xi = randn(N,n);
case 2
  Xi = unifrnd(-sqrt(3),sqrt(3),N,n);
otherwise
  error('Not implemented yet!!');
end
