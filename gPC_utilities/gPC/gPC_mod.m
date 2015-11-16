function um = gPC_mod(w,u,Psi)
%GPC_MOD returns a set of spectral gPC modal coefficients from a given realization ensemble
%associated to a discrete probability measure via some weights. The algorithm relies on 
%numerical quadrature
%   GPC_MOD(W,U,Psi), returns values of U in the modal space Um.
%   The first argument is the column vector of distribution weights.
%   The second argument are the realizations of the RV(s) U. Each row of U is a
%   realization and each column a RV.
%   The third argument represents the gPC polynomial basis evaluated at some points
%   (realizations of the underlying RV(s) seed(s)). Each row of Psi is a realization and
%   each column a polynomial member.
%
%   Remark: vector/matrice arguments must be properly ordered

[Nq,P] = size(Psi); [Nq,n]  = size(u);
um = zeros(P,n);

for i=1:n,
    for j=1:P,
        um(j,i) = Numquad(w,u(:,i).*Psi(:,j));
    end
end
