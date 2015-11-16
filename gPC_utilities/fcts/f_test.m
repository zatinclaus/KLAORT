function res=f_test(x,cas)

res=zeros(size(x));

if strcmp(cas,'bw_facing_step')==1
   shift=-0.4;
   level_low = 0; level_high = 1;
   res=step_func(x,shift,level_high,level_low);
end;

if strcmp(cas,'bw_facing_x2step')==1
   shift1=-0.4; shift2=0.6;
   level_low = 0; level_mid = 0.5; level_high = 1;
   res=step_func(x,shift1,level_mid,level_low)+step_func(x,shift2,level_mid,level_low);
end;

if strcmp(cas,'gas_dyn_like')==1
   shift1=-0.6; shift2=0;
   level_low = 0; level_mid = 1.0;
   res=step_func(x,shift1,level_low,level_mid)+step_func(x,shift2,level_low,level_mid).*(x-shift2).^2;
end;

if strcmp(cas,'sin')==1
   res=sin(pi*x);
end;

if strcmp(cas,'tan')==1
   res=tan(x);
end;

if strcmp(cas,'tanh')==1
   coeff = 6; shift = 0.2; x_shift = 0.1;
   res=tanh(coeff.*(x-x_shift))+shift;
end;

if strcmp(cas,'rat')==1
   delta = 1.4;
   res = 1./((x-delta).^2);
end;

if strcmp(cas,'LegMpoly')==1 % Legendre Monic polynomials
   index = 5;
   ab_leg = r_jacobi(index+1);
   res = orthopol(ab_leg,x,index);
end;

if strcmp(cas,'siampoly')==1 % Legendre Monic polynomials
   index1 = 0; index2=3; index3=10;
   ab_leg1 = r_jacobi(index1+1); ab_leg2 = r_jacobi(index2+1); ab_leg3 = r_jacobi(index3+1);
   res = orthopol(ab_leg1,x,index1)+ orthopol(ab_leg2,x,index2)+ 40.*orthopol(ab_leg3,x,index3);
end;

if strcmp(cas,'sin_dur')==1
   res=sin(2*pi*x);
end;

if strcmp(cas,'sin_aff')==1
   res=1+5*x+sin(2*pi*x);
end;

if strcmp(cas,'indic_sin')==1
   translation=0;
   res=indic_sin(x,translation,10);
end;
if strcmp(cas,'hat')==1
   a = -1./3; b = 1./3; level_low = 0; level_high = 1;
   res=step_func(x,a,level_low,level_high).*step_func(x,b,level_high,level_low);
end;

if strcmp(cas,'bihat')==1
   res=chapeau(x,-0.5,1)+chapeau(x,0.5,1);
end;

if strcmp(cas,'trihat')==1
   res=chapeau(x,-0.5,1)+chapeau(x,0.5,1)-chapeau(x,0.,1);
end;
if strcmp(cas,'cuvette')==1
   %une cuvette
   res=chapeau(x,-0.5,1)+chapeau(x,0.5,1)+chapeau(x,0.,1);
end;

if strcmp(cas,'cren')==1
   %un creneau special
   res=heavy(x,-0.5,0.);
end;
if strcmp(cas,'gauss')==1
   %y est gaussien => u est uniforme sur 0,1
   %u=1./2+1./2*erf(y*sqrt(2)/2);
   %y est gaussien => u est uniforme sur -1,1
   %u=erf(y*sqrt(2)/2);
   
   %Attention, on ne peut pas utiliser Clenshaw Curtis car pas def en -1
   %et 1
   res=erfinv(x)*2./sqrt(2);
end;
if strcmp(cas,'indic_aff')==1
   translation=-0.5;
   res=indic_aff(x,translation);
end;
if strcmp(cas,'indic_affp')==1
   translation=0.5;
   res=indic_affp(x,translation);
end;

if strcmp(cas,'exp1')==1
   res=exp(-x);
end;

if strcmp(cas,'exp10')==1
   %cas de gerritsma
   res=exp(-x*10);
end;
if strcmp(cas,'exp15')==1
   %cas de gerritsma
   res=exp(-x*15);
end;
if strcmp(cas,'exp100')==1
   %cas de gerritsma
   res=exp(-x*100);
end;

if strcmp(cas,'rampe')==1
   %cas rampe (pr??vision Burgers)
   res=rampe(x,0.25);
end;

if strcmp(cas,'expPl')==1
   res = exp(-200*x.^2)-x;
end;

if strcmp(cas,'lineshift')==1
   shift = 1000;
   res = x+shift;
end;
