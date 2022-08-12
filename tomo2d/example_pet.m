
%
fin_path = '~/Work/fastfins.m';
warning('off')
rmpath(genpath(fin_path));
warning('on')
addpath(genpath([fin_path '/fastfins/']));
addpath(genpath([fin_path '/tomo2d']));

fastfins_check_solvers()

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else

    tomo_opt = tomo_options('mesh_size', 64, 'angle', pi*2, 'n_sourc', 100, 'Is',  50);

    %inv_opts = inverse_options('cov_type', 'MRF', 'mean', 0, 'k', 100, 'cond', [1, 1, 0], 'sigma', 4, ...
    %    'test_type', 'CF', 'test_base', 2, 'test_range', -4);

    inv_opts = inverse_options('cov_type', 'GP', 'mean', 0, 'scale', 200, 'power', 1, 'sigma', 2, ...
        'test_type', 'CF', 'test_base', 2, 'test_range', -5);

    [model, obs, prior] = setup_tomo(tomo_opt, inv_opts, '');

    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = zeros(prior.dof, 1);
    vmap = get_map_2020a(mlpost, mvhess, vini);
    umap = matvec_prior_L(prior, vmap) + prior.mean_u;

    figure
    subplot(1,2,1)
    surf(reshape(prior.trueu, tomo_opt.mesh_size, tomo_opt.mesh_size), 'edgecolor', 'none')
    subplot(1,2,2)
    surf(reshape(umap, tomo_opt.mesh_size, tomo_opt.mesh_size), 'edgecolor', 'none')
    
end

return

mcmc_opt = mcmc_options('nstep', 1E2, 'sbatch', 1, 'lis_flag', true, 'sigma', -1);
hess = @(v) eigen_PPFisher(model, obs, prior, v, 1E-2, obs.n_data);
out_hmala = hmala(mlpost, hess, vmap, mcmc_opt);
Hg = out_hmala.cross_g/out_hmala.nsample;

%{
ind = 1:1E2:out_hmala.nsample;

Hvar = zeros(prior.dof);
Hmean = zeros(prior.dof);
for i = 1:numel(ind)
    v = out_hmala.samples(:,i);
    u = matvec_prior_L(prior, v) + prior.mean_u;
    sol = forward_solve(model, u);
    Hv = hessian(model, prior, sol, sol.d-obs.data);
    Hmean = Hmean + Hv;
    Hvar = Hvar + Hv*Hv;
end

Hmean = Hmean/numel(ind);
Hvar = Hvar/numel(ind) - Hmean.^2;

subplot(2,2,1)
surf(Hg)

subplot(2,2,2)
surf(Hmean)

subplot(2,2,3)
surf(Hvar)
%}

%{
v = vmap;
u = matvec_prior_L(prior, v) + prior.mean_u;
sol = forward_solve(model, u);

tol = 1E-3;
gp = zeros(prior.dof);
gm = zeros(prior.dof);
for i = 1:numel(v)
    vp = v;
    vp(i) = vp(i) + tol;
    vm = v;
    vm(i) = vm(i) - tol;
    [~,~,~,gp(:,i)] = mlpost(vp);
    [~,~,~,gm(:,i)] = mlpost(vm);
end
Hd = (gp - gm)/(2*tol);

Hv = hessian(model, prior, sol, sol.d-obs.data);
%}


