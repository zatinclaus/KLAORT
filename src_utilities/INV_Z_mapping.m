function [zm] = INV_Z_mapping(z,zmin,zmax)
% Linear mapping from [zmin,zmax] -> [-1,1]
zm = (2.*z - (zmax+zmin))./(zmax-zmin);
end