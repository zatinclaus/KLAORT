function [Ys] = CSSmat(X,Y,XX,dim,varargin)
%
%  CSSmat(X,Y,XX,dim,P) Cubic Smoothing Spline for matrices Y.
%  uses the csaps cubic smoothing spline routine.
%  This function returns the smoothed data matrix
%  The smoothing is done along the matrix dimension dim.
%  The original grid is X and the new grid is XX (can be the same).
%  P controls the smoothing and must be between 0 (linear fit)
%  and 1 (normal spline fit), otherwise it is self assigned.
%  [Ys] = CSSmat(X,Y,XX,dim,P);

if nargin == 4,
  fprintf(1,'The amount of smoothing P will be automatically chosen to be P=0.9\n');
  P = 0.9;
else
 P = varargin{:};
end

perm = [dim:max(ndims(Y),dim) 1:dim-1];
Y = permute(Y,perm);

[l,c] = size(Y);
Ys = zeros(l,length(XX));

fprintf(1,'The smoothing will be done along the dimension with n=%d components.\n',c);

for i=1:l,
      Ys(i,:) = Csaps(X,Y(i,:),P,XX);
end
