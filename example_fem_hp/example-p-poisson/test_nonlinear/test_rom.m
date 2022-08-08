

renew_setup = false;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    p_rate = 3;
    epsilon = 1E-6;
    tol = 1E-4;
    bnd_funcs = {@(p) find(abs(p(:,2))<tol)};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1E3;
vs = randn(prior.dof, n);
%
rom_opts.deim_reg_factor = 1.2;
rom_opts.deim_state_tol  = 1E-4;
rom_opts.pod_state_tol   = 1E-5;
%
states = zeros(model.mesh.dof,n);
tic;
for i = 1:n
    [~,~,~,~,sol] = minus_log_post(model, obs, prior, vs(:,i));
    states(:,i) = sol.state;
end
toc
weights = ones(1,n)/n;
%
prior_KL = basis_KL(prior, 1-1E-2);
%
sub_vs = prior_KL.chol2w'*vs;
%
tic;
rom = setup_p_poisson_rom(model, prior_KL, sub_vs, states, weights, rom_opts);
toc
rom.res_tol = 1E-5;

tic;
for i = 1:20
    sol = forward_solve(model, prior_KL.basis*sub_vs(:,i)+prior_KL.mean_u);
    sta_f(:,i) = sol.state;
end
toc
tic;
for i = 1:20
    %
    solr = rom_solve(rom, sub_vs(:,i));
    sta_r(:,i) = solr.state;
end
toc
tic;
for i = 1:20
    %
    sols = forward_solve_redu(model, rom.states, prior_KL.basis*sub_vs(:,i)+prior_KL.mean_u);
    sta_s(:,i) = sols.state;
end
toc
