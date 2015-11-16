function [ y ] = sinpuls(t,Tw)
%SINPULS Sampled aperiodic half-sine generator.
%   SINPULS(T) generates samples of a continuous, aperiodic,
%   unity-height half-sine at the points specified in array T, centered
%   about T=0.  Note that the
%   interval of non-zero amplitude is defined to be open on the right,
%   i.e., SINPULS(-0.5)=1 while SINPULS(0.5)=0.
%
%   SINPULS(T,W) generates a half-sine of width W.

% error(nargchk(1,2,nargin,'struct'));
% if nargin<2, Tw=1;   end

y = abs(t)<Tw/2-eps;
y = y.*sin(pi.*(t+(Tw./2))./Tw);
