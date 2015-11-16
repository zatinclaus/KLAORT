function [ym,ind] = findpeak(x,y)
% find local maxima and minima and return them along 
% with the indices where they occur
cpzr = 1.0e-09 ;  % computational zero
nx = max(size(x,1),size(x,2)) ;
ny = max(size(y,1),size(y,2)) ;
if (nx~=ny)
    disp('vectors x, y should be of equal size') ;
    return ;
else
    nn = nx ;
end
% sort x-values 
[xst, indst] = sort(x) ;
yst = y(indst) ;
% eliminate possible duplicates
xsr = zeros(nn,1) ;
ysr = zeros(nn,1) ;
indr = zeros(nn,1) ;
k = 1 ;
xsr(1) = xst(1) ;
ysr(1) = yst(1) ;
for i=2:nn
   if (abs(xst(i)-xst(i-1))<cpzr)
      continue ;
   else
      k = k + 1 ;
      xsr(k) = xst(i) ;
      ysr(k) = yst(i) ;
      indr(k) = indst(i) ;
   end    
end 
n = k ;
xs = xsr(1:k) ;
ys = ysr(1:k) ;
inds = indr(1:k) ;
% -----------------------
slpback = NaN*ones(n,1) ;
slpfore = NaN*ones(n,1) ;
for i=2:n
   slpback(i) = (ys(i)-ys(i-1))/(xs(i)-xs(i-1)) ;    
end
for i=1:n-1
   slpfore(i) = (ys(i+1)-ys(i))/(xs(i+1)-xs(i)) ;  
end

inda = zeros(n,1) ;
ya = zeros(n,1) ;
k = 0 ;
for i=2:n-1
   if (slpback(i)*slpfore(i) <= 0)
       k = k + 1 ;
       inda(k) = inds(i) ;
       ya(k) = y(inds(i)) ;
   end    
end     
ym = ya(1:k) ;
ind = inda(1:k) ;
