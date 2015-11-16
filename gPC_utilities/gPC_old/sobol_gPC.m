function sobol = sobol_gPC(dim,order,multindex,x,varargin)
%SOBOL_GPC Sobol coefficients matrix up to second-order cross terms
%   obtained from gPC modal decomposition
%   SOBOL_GPC(X), if X is a vector, returns a vector.  For matrices,
%   where each row is a gPC modal coefficient, and each column a variable,
%   SOBOL_GPC(X) returns a matrix of coefficients.
%   SOBOL_GPC(X,P), where P is a vertical vector of gPC basis norms.
%   Default value for P is the identity vector (useful when gPC approximation
%   basis is built as orthonormal polynomials).
%
%   Remark: the mean (i.e. first gPC modal coefficient) is removed from each
%   column before calculating the result.

nin = nargin;

if nin==0
    error('MATLAB:sobol_gPC:NotEnoughInputs','Not enough input arguments.');
end
if nin>6
    error('MATLAB:sobol_gPC:TooManyInputs', 'Too many input arguments.');
end
if ~ismatrix(x)
    error('MATLAB:sobol_gPC:InputDim', 'Inputs for x must be at max 2-D.');
end

scalarxy = false; % cov(scalar) is an ambiguous case
[m,n] = size(x);

if nin == 4, % only parameters and x are provided
    poly_norm = ones(m,1); % default polynomial norm taken to 1.0
elseif nin == 5, % table of coefficients for S_ij is provided
    poly_norm = ones(m,1); % default polynomial norm taken to 1.0
    tab = varargin{1};
else % i.e. nin = 6
    poly_norm = varargin{2};
    if ~ismatrix(poly_norm)
        error('MATLAB:cov_gPC:InputDim', 'Second input, i.e. for polynomial norm must be 1-D.');
    elseif length(poly_norm(1,:))>1,
        error('MATLAB:cov_gPC:InputDim', 'Second input, i.e. for polynomial norm must be 1-D vertical vector.');
    end
    if length(x(:,1)) ~= length(poly_norm),
        error('MATLAB:cov_gPC:XYlengthMismatch', 'The number of rows in x and polynomial norm must match.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isempty(x);
%    if (m==0 && n==0)
%       xy = NaN(class(x));
%    else
%       xy = NaN(n,class(x));
%    end
%    return;
% end
%
% if m == 1  % One observation
%    % For single data, unbiased estimate of the covariance matrix is not defined.
%    % Return the second moment matrix of the observations about their mean.
%    xy = zeros(n,class(x));
% else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COV = cov_gPC(x,poly_norm);
VAR = diag(squeeze(COV))'; clear COV;

x_minus_mean = x(2:end,:); x = x_minus_mean; clear x_minus_mean;
multindex = squeeze(multindex(2:end,:));

poly_norm = poly_norm(2:end);

if order>1,
    n_sob = PCnumbterms(dim,2) - PCnumbterms(dim,1);
else
    n_sob = PCnumbterms(dim,1)-1;
end

sobol = zeros(n_sob,n);

for i=1:PCnumbterms(dim,1)-1, % contributions from non-coupled linear terms
    sobol(i,:) = x(i,:).^2.*poly_norm(i);
end

if order>1,
    for i=PCnumbterms(dim,1):PCnumbterms(dim,order)-1, % contributions from non-coupled non-linear terms
        ind_tmp = find(multindex(i,:));
        if length(ind_tmp) == 1,
            sobol(ind_tmp,:) = sobol(ind_tmp,:) + x(i,:).^2.*poly_norm(i);
        end
    end
    clear ind_tmp;
    
    for i=PCnumbterms(dim,1):PCnumbterms(dim,order)-1, % contributions from second-order coupled non-linear terms
        
        ind_tmp = find(multindex(i,:));
        if length(ind_tmp) == 2,
            j = find(tab(:,1) == ind_tmp(1) & tab(:,2) == ind_tmp(2));
            if isempty(j), warning('sobol_gPC.m: problem with Sobol indices ordering.'); end
            sobol(tab(j,3),:) = sobol(tab(j,3),:) + x(i,:).^2.*poly_norm(i);
        end
    end
end

sobol = sobol./(repmat(VAR,n_sob,1));

end