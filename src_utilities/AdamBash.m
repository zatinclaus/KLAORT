function [U_new] = AdamBash(Delta_t,U_old,Fn,Fnm1)

U_new = U_old + Delta_t.*(1.5.*Fn -0.5 .*Fnm1);
 