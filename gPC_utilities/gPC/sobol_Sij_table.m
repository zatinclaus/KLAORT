function [table] = sobol_Sij_table(dim,n)
% this function creates a table correspondance that facilitates the sorting
% and putting away of second-order cross Sobol' indices

table = zeros(n-dim,3);
count = 0;
for i=1:dim,
    for j=i+1:dim,
        count = count+1;
        table(count,1:2) = [i j];
        table(count,3) = dim+count;
    end
end

end

