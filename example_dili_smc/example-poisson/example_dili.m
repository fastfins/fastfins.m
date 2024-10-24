

%% define the target density 
ml_target = @(v) minus_log_post(model, obs, prior, v);
% define the matrix vector product with the Hessian or the Gauss-Newton
% Hessian
hessian = @(sol,tol,nmax) eigen_PPFisher(model, obs, prior, sol, tol, nmax);


% find map
mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
vini = randn(prior.dof, 1);
vmap = get_map_2023a(ml_target, mvhess, vini);

% options for building the LIS
% type help build_lis for details (not written yet)
% here the options are using adaptive MCMC ('method', 'DILI')
% using eigen decomposition ('hess_type', 'Eig')
% maximum number of parameter pints for evaluating the Hessian ('max_hess', 100)
% maximum rank for each local Hessian ('local_max_dim', 50)
% the storage rank for saving the global Hessian ('global_max_dim', 200)
% See the DILI paper for details of the iterative construction

%%
lis_opt_dili = lis_options('method', 'DILI', 'hess_type', 'Eig', ...
    'max_hess', 100, 'min_hess', 10, 'local_max_dim', 50, 'global_max_dim', 200);

% build LIS
[P_dili,S_dili,gsvd_dili,out_dili] = build_lis(ml_target, hessian, vmap, lis_opt_dili);
stat = out_dili.stat;

lis_dili.V = P_dili;
lis_dili.s = S_dili;
lis_dili.r = size(P_dili,2);

%% run DILI MCMC using the LIS
opt = mcmc_options('proposal', 'MALA', 'nstep', 1E4, 'sbatch', 1, 'sigma', -3);
tic
out_dili = dili(ml_target, P_dili, stat, vmap, opt);
toc

%%
% plots and examples of how to extract results from out_dili
figure
plot(out_dili.mllkd)
title('minus log likelihood')

figure
plot(out_dili.mlpt)
title('minus log posterior')

figure
plot(out_dili.sub_samples(1,:))
title('the first element of subspace markov chain')

figure
plot(out_dili.samples(1,:))
title('the first element of full space markov chain')