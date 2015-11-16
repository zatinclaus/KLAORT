function [ ] = KDE2D_plot(data,fig_number,b_and_w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[dim,N] = size(data); dim_sq = dim^2;

p = kde(data,'rot');
mean(p)
count = 0;

figure(fig_number); hold on; gray_dark = 0.7.*ones(3,1);

for i=1:dim,
    for j=1:dim,
        count = count+1;
        if (j<i)
            subplot(dim,dim,count);
            pm = marginal(p,[i;j]);
            if b_and_w,
                [H,X,Y] = hist(pm,100); contour(X,Y,H,8,'lines','-','color',gray_dark,'linewidth',0.5); axis image; axis square% colorbar;
            else
                [H,X,Y] = hist(pm,100); contour(X,Y,H,10,'lines','-','linewidth',1.5); axis image; axis square% colorbar;
            end
            hold on;
        end
        if (j==i),
            subplot(dim,dim,count);
            if b_and_w,
                pm = marginal(p,[i]); plot(pm,'-k'); axis square; hold on;
            else
                pm = marginal(p,[i]); plot(pm,'-r'); axis square; hold on;
            end
        end
    end
end

