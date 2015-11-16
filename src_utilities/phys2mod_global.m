function [U_mod] = phys2mod_global(Phi_nodes_weighted, U_phys, M_norm, nq)

U_mod = [Phi_nodes_weighted*U_phys(1:nq,:);Phi_nodes_weighted*U_phys(nq+1:2*nq,:)]; % faster (about /2) than concatenating with cat
U_mod = U_mod./M_norm;
