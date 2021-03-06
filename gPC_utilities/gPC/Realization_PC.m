function [Uf] = Realization_PC(Um,Psi,N,T,P)
% This function computes N individual time T (or space X) random
% realizations given the first P random GPC coefficients Um(t/x) and P
% polynomials Psi(xi) of the process. size(Psi) = [N P] and size(Um) = [T P]
% and size(Uf) = [N T]

% [r,c]  = size(Psi); [rr,cc] = size(Um);
% 
% if (r~=N || c~=P)
%     error('The gPC polynomial basis must have dimension [N,P]');
% end
% 
% if (rr~=T || c~=P)
%     error('The gPC modal coefficients must have dimension [T,P]');
% end


Um  = repmat(Um(1:T,1:P),[1 1 N]); Psi = repmat(Psi(1:N,1:P),[1 1 T]);

Uf = sum((permute(Um,[3 1 2])).*(permute(Psi,[1 3 2])),3);