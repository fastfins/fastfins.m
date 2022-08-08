load_dir;

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    p_rate = 3;
    epsilon = 1E-6;
    tol = 1E-4;
    bnd_funcs = {@(p) find(abs(p(:,2))<tol); @(p) find(abs(p(:,1))<tol)};
    %
    obs_locs    = [linspace(0.1,0.9,20)', 0.05*ones(20,1)];
    force_func  = @(x1, x2) 1;
    %
    model_opts = hp_options('h', 1/50, 'poly_order', 2, 'quad_order', 10, ...
        'xyratio', 0.2, 'obs_locs', obs_locs, 'bnd_funcs', bnd_funcs, 'force_func', force_func, ...
        'nl_solver', 'line_search', 'res_tol', 1E-10, 'exp_param', true, 'nl_disp', 0);
    %
    inv_opts = inverse_options('s2n', 50, 'cov_type', 'GP', 'mean', 0, ...
        'power', 1, 'scale', 0.5, 'sigma', 1);
    %
    [model, obs, prior] = setup_p_poisson(model_opts, inv_opts, p_rate, epsilon, []);
    
    solt = forward_solve(model, obs.true_u);
    figure
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    figure
    plot(model.mesh.nodes(model.pa_bnd_mesh.node_indices, 1), obs.true_u)
    
    
    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = randn(prior.dof, 1);
    vmap = get_map_2020a(mlpost, mvhess, vini);
    umap = matvec_prior_L(prior, vmap) + prior.mean_u;
    solm = forward_solve(model, umap);
    
    figure
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solm.state, 'edgecolor', 'none')
    figure
    plot(model.mesh.nodes(model.pa_bnd_mesh.node_indices, 1), umap)
    
    figure
    plot(obs.data, 'o')
    hold on
    plot(solm.d, 'x')
    plot(solt.d, '-')
    
end

% test the qoi
tic;
for i = 1:100
    u = matvec_prior_L(prior, randn(prior.dof,1)) + prior.mean_u;
    sol = forward_solve(model, u);
    qois(i) = sol.qoi;
end
toc



