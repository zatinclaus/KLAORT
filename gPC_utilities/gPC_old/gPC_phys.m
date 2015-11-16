function [ux] = gPC_phys(um,Psi)
%GPC_PHYS returns a set of realizations of a RV from its spectral gPC modal decomposition
%   GPC_PHYS(Um,Psi), returns values of U in the space of realizations Ux.
%   The first argument are the gPC modal coefficients. It may be a line vector or a matrix.
%   In this case, each row is a gPC modal coefficient, and each column a different RV.
%   The second argument represents the gPC polynomial basis evaluated at some points
%   (realizations of the underlying RV(s) seed(s)). Each row of Psi is a
%   realization and each column a polynomial member.
%
%   Remark: vector/matrice arguments must be properly ordered


% This function computes Nq individual time T (or space or number of RVs) random
% realizations given the first P random GPC coefficients um(t/x) and P
% polynomials Psi(xi) of the process. size(Psi) = [N P] and size(um) = [T P]
% and size(ux) = [N T]

[Nq,P]  = size(Psi); [T,P] = size(um);

um  = repmat(um(1:T,1:P),[1 1 Nq]);
Psi = repmat(Psi(1:Nq,1:P),[1 1 T]);
um  = permute(um,[3 1 2]);
Psi = permute(Psi,[1 3 2]);

ux = sum(um.*Psi,3);