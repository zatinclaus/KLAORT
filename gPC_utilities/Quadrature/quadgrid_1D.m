function [x,w,N] = quadgrid_1D(quad_type,Nq,varargin)
%QUADGRID_1D quadrature points and weights.
%
% [x,w,N] = quadgrid_1D(quad_type,Nq)
% [x,w,N] = quadgrid_1D(quad_type,Nq,'nested')
%
% The routine computes the location of the quadrature points and the
% corresponding quadrature weights of Clenshaw-Curtis, Fejer or
% Gauss-Legendre quadratures in 1D in [-1,1]
%
% INPUTS
%    quad_type     type of quadrature (CC, FJ, GL, GH, EQ or MC)
%    Nq            number of wished quadrature points
%
% OPTIONAL INPUTS
%    'nested'      adjusts N to get nested grids (only valid for CC and FJ)
%
% OUTPUTS
%    x     N quadrature points
%    w     N quadrature weights with sum(w) = 1
%    N     number of quadrature points retained


% default values
nested=[];

for i=1:1:length(varargin)
    switch lower(varargin{i})
        case 'nested'
            nested = 1;
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if strcmp(quad_type,'CC')==1 % Clenshaw-Curtis type quadrature
    
    if nested,
        t=[1 2.^((2:1:13)-1)+1];
        [dum1,t_index] = min(abs(t-Nq));
        N = t(t_index);
    else
        N = Nq;
    end
    
    [x,w] = clencurt(N-1); x = sort(x); w = (w./2)';
    
    
elseif strcmp(quad_type,'GL')==1 % Gauss-Legendre type quadrature
    
    N = Nq;
    [x,w] = GLNodeWt(N); w = w./2;
    
elseif strcmp(quad_type,'GH')==1 % Gauss-Hermite type quadrature
    
    N = Nq; ab = r_hermite(N);
    [xw] = gauss(N,ab); x = xw(:,1).*sqrt(2); w = xw(:,2)./sqrt(pi);
    
elseif strcmp(quad_type,'FJ')==1  % Fejer type quadrature
    
    if nested,
        t=2.^(1:1:12)-1; [dum1,t_index] = min(abs(t-Nq));
        N = t(t_index);
    else
        N = Nq;
    end
    
    xw = fejer(N);
    x = xw(:,1); w = xw(:,2)./2; clear xw;
    
elseif strcmp(quad_type,'EQ')==1  % Equidistant type quadrature
    
    N = Nq;
    x = linspace(-1,1,N)';
    w = (1./N).*ones(size(x));
    
elseif strcmp(quad_type,'MC')==1  % MC uniformly distributed type quadrature
    
    N = Nq;
    x = sort(unifrnd(-1,1,N,1));
    w = (1./N).*ones(size(x));
    
else error('Quadrature not implemented yet!');
    
end

if nested,
    disp(['Number of quadrature points retained for nested grid Nq=',num2str(N)]);
end