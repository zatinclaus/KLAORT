function [U_phys] = mod2phys_global(Phi_nodes,U_mod,P)

U_phys = [Phi_nodes*U_mod(1:P,:);Phi_nodes*U_mod(P+1:2*P,:)]; % faster (about /2) than concatenating with cat