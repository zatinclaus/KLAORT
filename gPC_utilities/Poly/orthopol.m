% ORTHOPOL  Orthogonal polynomial at location x, based on the first
% n three-term recurrence coefficients.
%
%   fx=ORTHOPOL(ab,x,n) evaluates the orthogonal MONIC polynomial of order n+1 at the location x.
%   The alpha-coefficients must be stored in the first column and the beta-coefficients
%   in the second column, of the nx2 array ab.
%   The call fx=ORTHOPOL(ab,x) is the same as fx=ORTHOPOL(ab,x,length(ab(:,1))).
%
function fx = orthopol(ab,x,n)
if nargin<3, n=length(ab(:,1)); end
if(n<0) error('parameter(s) out of range'), end

fx = 0.0; fxm1 = 1.0; fxm2 = 0.0;

if ~n,
  fx = ones(size(x));
else
  for k=1:n,
    fx = ((x-ab(k,1)).*fxm1 - ab(k,2).*fxm2); fxm2 = fxm1; fxm1 = fx;
  end
end
