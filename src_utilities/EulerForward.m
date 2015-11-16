function [U_new] = EulerForward(Delta_t,U_old,Fn)

U_new = U_old + Delta_t.*Fn;
 