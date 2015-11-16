function xy = cov_gPC(x,varargin)
%COV_GPC Covariance matrix obtained from gPC modal decomposition
%   COV_GPC(X), if X is a vector, returns the variance.  For matrices,
%   where each row is a gPC modal coefficient, and each column a variable,
%   COV_GPC(X) is the covariance matrix.  DIAG(COV_GPC(X)) is a vector of
%   variances for each column, and SQRT(DIAG(COV_GPC(X))) is a vector
%   of standard deviations. COV_GPC(X,P), where P is a vertical vector of gPC
%   basis norms. Default value for P is the identity vector (useful when
%   gPC approximation basis is built as orthonormal polynomials).
%
%   Remark: the mean (i.e. first gPC modal coefficient) is removed from each
%   column before calculating the result.

nin = nargin;

if nin==0
   error('MATLAB:cov_gPC:NotEnoughInputs','Not enough input arguments.');
end
if nin>2
   error('MATLAB:cov_gPC:TooManyInputs', 'Too many input arguments.');
end
if ~ismatrix(x)
   error('MATLAB:cov_gPC:InputDim', 'Inputs for x must be at max 2-D.');
end

scalarxy = false; % cov(scalar) is an ambiguous case
[m,n] = size(x);

if nin == 1, % only x is provided 
    poly_norm = ones(m,1); % default polynomial norm taken to 1.0
else % i.e. nin = 2
   poly_norm = varargin{1};
   if ~ismatrix(poly_norm)
      error('MATLAB:cov_gPC:InputDim', 'Second input, i.e. for polynomial norm must be 1-D.');
   elseif length(poly_norm(1,:))>1,
      error('MATLAB:cov_gPC:InputDim', 'Second input, i.e. for polynomial norm must be 1-D vertical vector.');
   end
   if length(x(:,1)) ~= length(poly_norm),
      error('MATLAB:cov_gPC:XYlengthMismatch', 'The number of rows in x and polynomial norm must match.');
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(x);
   if (m==0 && n==0)
      xy = NaN(class(x));
   else
      xy = NaN(n,class(x));
   end
   return;
end

if m == 1  % One observation   
   % For single data, unbiased estimate of the covariance matrix is not defined.
   % Return the second moment matrix of the observations about their mean.
   xy = zeros(n,class(x));
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_minus_mean = x(2:end,:); x = x_minus_mean; clear x_minus_mean;
poly_norm = poly_norm(2:end);

xy = x'*(x.*repmat(poly_norm,1,n));
   
end
