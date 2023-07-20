% This script demonstrates running the HMALA algorithm, using the inverse
% posterior Hessian at the MAP to precondition the Langevin equation
% Tiangang Cui, 17/May/2014

ml_target = @(v) minus_log_post(model, obs, prior, v);
matvec_hessian = @(HI, dv) matvec_PPFisher(model, obs, prior, HI, dv);
hessian   = @(HI) eigen_PPFisher(model, obs, prior, HI, 1E-4, obs.n_data-1);

% find map
%vmap1 = get_map_v1(model, obs, prior, zeros(prior.DoF, 1));
%vmap2 = get_map_v2(ml_target, matvec_hessian, zeros(prior.DoF, 1));


opt = mcmc_options('proposal', 'MALA', 'nstep', 1E4, 'sbatch', 2, 'sigma', -1.5, ...
    'lis_flag', false);

tic
out_hmala = hmala(ml_target, hessian, vmap, opt);
toc

