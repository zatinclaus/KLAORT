function [f] = f_ramp(time,ti,tf)
  fmin = 0; fmax = 1;
  f = ((fmax-fmin)./(tf-ti)).*(time-ti) + fmin;