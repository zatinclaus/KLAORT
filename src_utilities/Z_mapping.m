function [zm] = Z_mapping(z,zmin,zmax)
% Linear mapping from  [-1,1] -> [zmin,zmax]
zm = (zmax-zmin).*z./2 + (zmax+zmin)./2;