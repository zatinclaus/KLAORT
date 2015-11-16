function [um,gPC_mse] = gPC_mod_PSP(w,u,Psi)
%GPC_MOD_PSP returns a set of spectral gPC modal coefficients from a given realization cubature ensemble
%associated to a discrete probability measure via some weights. The pseudo-spectral projection algorithm
%relies on numerical quadrature. It also returns an estimate of the std of the projection error term
%   GPC_MOD_PSP(W,U,Psi), returns values of U in the modal space Um and
%   mean square error of the gPC approximation on the set of data points
%   The first argument is the column vector of distribution weights.
%   The second argument are the realizations of the RV(s) U. Each row of U is a
%   realization and each column a RV.
%   The third argument represents the gPC polynomial basis evaluated at some points
%   (realizations of the underlying RV(s) seed(s)). Each row of Psi is a realization and
%   each column a polynomial member.
%
%   Remark: vector/matrice arguments must be properly ordered

[~,P] = size(Psi); [Nq,n]  = size(u); um = zeros(P,n);

for i=1:n,
    for j=1:P,
        um(j,i) = Numquad(w,u(:,i).*Psi(:,j));
    end
end
UM  = repmat(um(1:P,1:n)',[1 1 Nq]); PSI = repmat(Psi(1:Nq,1:P),[1 1 n]);
ugPC = sum((permute(UM,[3 1 2])).*(permute(PSI,[1 3 2])),3);

gPC_mse = squeeze(sum((ugPC-u).^2,1));
% [l,~] = size(gPC_mse); if l==1, gPC_mse = gPC_mse'; end

end
