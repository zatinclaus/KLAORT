function [um, gPC_mse] = gPC_mod_LASSO(u,Psi,CODE,varargin)
%GPC_MOD_LASSO returns a set of spectral gPC modal coefficients from a given realization ensemble
%associated to a continous probability measure via a Regularized least-squares regression using lasso.m
%It also returns an estimate of the std of the regression error term.
%GPC_MOD_LASSO(U,PSI,CODE,VARARGIN), returns values of U in the modal space Um,
%and mean square error: gPC_mse
%   The first argument are the realizations of the RV(s) U. Each row of U is a
%   realization and each column a different RV.
%   The second argument represents the gPC polynomial basis evaluated at some points
%   (realizations of the underlying RV(s) seed(s)). Each row of Psi is a realization and
%   each column a polynomial member.
%   CODE is used mainly to determine the parallelism flag
%   Remark: vector/matrice arguments must be properly ordered
weight_flag = false;
if nargin == 4,
    weight = varargin{:}; weight_flag = true;
end

verbose = 1; viewvose = 0;
[Nq,P] = size(Psi); [Nq,n]  = size(u);
m = length(CODE.n_QoI); um = zeros(P,n); % and not zeros(P,m);

CV_fold_nbr = min(Nq,10); % Nq; % number of folds for cross-validation
epsilon_tol = 1e-4; % 1e-4 by default (lasso cpu time is sensitive to this parameter)
MC_Reps_nbr = 1;
lamb_rat = 1e-4;

% If one wants to impose predefined Lambda vector values:
% Lambda_vec = logspace(-5,-1,100);
% [b,fitinfo] = lasso(...,'Lambda',Lambda_vec); % to be used as is

if CODE.parallel,
    opts = statset('UseParallel',true);
end

if verbose,
    disp(['gPC_mod_LASSO: for N=',num2str(Nq),', samples:']);
end
for i=1:m,
    ii = CODE.n_QoI(i);
    if weight_flag,
        tic;
        if CODE.parallel,
            [b,fitinfo] = ...
                lasso(Psi,u(:,ii),'CV',CV_fold_nbr,...
                'Options',opts,'RelTol',epsilon_tol,'MCReps',MC_Reps_nbr,'Weights',weight,'LambdaRatio',lamb_rat);
        else
            [b,fitinfo] = ...
                lasso(Psi,u(:,ii),'CV',CV_fold_nbr,...
                'RelTol',epsilon_tol,'MCReps',MC_Reps_nbr,'Weights',weight,'LambdaRatio',lamb_rat);
        end
        time = toc;
    else
        tic;
        if CODE.parallel,
            [b,fitinfo] = ...
                lasso(Psi,u(:,ii),'CV',CV_fold_nbr,...
                'Options',opts,'RelTol',epsilon_tol,'MCReps',MC_Reps_nbr,'LambdaRatio',lamb_rat);
        else
            [b,fitinfo] = ...
                lasso(Psi,u(:,ii),'CV',CV_fold_nbr,...
                'RelTol',epsilon_tol,'MCReps',MC_Reps_nbr,'LambdaRatio',lamb_rat);
        end
        time = toc;
    end
    
    % ind = fitinfo.Index1SE; % for a compromise with value retained within
    %    1 std error of the minimum
    ind = fitinfo.IndexMinMSE; % for the absolute minimum
    
    if verbose,
        disp(['Retained value of lambda: l=',num2str(fitinfo.Lambda(ind)),...
            ', for non-zero DoF: DF=',num2str(fitinfo.DF(ind)),...
            ', for mean squ. error: MSE=',num2str(fitinfo.MSE(ind)),', in: ',num2str(time),'(s)']);
        if viewvose,
            lassoPlot(b,fitinfo,'PlotType','Lambda','XScale','log');
            lassoPlot(b,fitinfo,'PlotType','CV');
            set(gca,'YScale','log'); % Use a log scale for MSE to see small MSE values better
            pause(0.1);
        end
    end
    
    um(1,ii) = fitinfo.Intercept(ind);
    um(2:end,ii) = b(2:end,ind);
    gPC_mse = fitinfo.MSE(ind);
    %     size(fitinfo.Lambda)
    %     [min(fitinfo.Lambda) max(fitinfo.Lambda)]
    %     figure(77); semilogx(fitinfo.Lambda,ones(size(fitinfo.Lambda)),'.','markersize',10,'color',rand(3,1)); hold on;
    
end

end
