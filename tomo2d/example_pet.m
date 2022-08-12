
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

    tomo_opt = tomo_options('mesh_size', 32, 'angle', pi*2, 'n_sourc', 50, 'n_detec', 100, 'd_width', 1, 'Is',  100);

    %inv_opts = inverse_options('cov_type', 'MRF', 'mean', 0, 'k', 5, 'cond', [1, 1, 0], 'sigma', 2, ...
    %    'test_type', 'CF', 'test_base', 2, 'test_range', -7);

    inv_opts = inverse_options('cov_type', 'GP', 'mean', 0, 'scale', 50, 'power', 1, 'sigma', 1, ...
        'test_type', 'CF', 'test_base', 2, 'test_range', -7);

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

