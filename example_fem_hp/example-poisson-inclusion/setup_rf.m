load_dir;

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    
    test_case = 'scalar';
    sq_flag = false;
    
    % boundary conditions
    tol = 1E-4;
    bnd_funcs = {@(p) find(abs(p(:,2))<tol), @(p) find(abs(p(:,2)-1)<tol)};
    bc_types  = {'essential', 'essential'};
    bc_funcs  = {@(x1,x2)-(cos(x1.*pi.*2)-1)/sqrt(2), @(x1,x2)cos(x1.*pi.*2)-1};
    
    % flux function for qoi
    flux_fx = @(x) uniform_fun (x, 0, 1);
    flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
    flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);
    
    % obs locations
    gx  = linspace(1/8, 7/8, 7);
    gy  = linspace(1/8, 7/8, 7);
    [xx0,yy0] = meshgrid(gx, gy);
    xx0 = xx0(:);
    yy0 = yy0(:);
    
    obs_locs  = unique([xx0, yy0], 'rows');
    
    %
    model_opts = hp_options('h', 1/32, 'poly_order', 1, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);
    
    inv_opts = inverse_options('s2n', 20, 'cov_type', 'MRF', 'mean', 0, ...
        'gamma', 20, 'cond', [0.2, 1, -1], 'sigma', 1);
    
    [model, obs, prior] = setup_poisson(model_opts, inv_opts, test_case, []);
    
    solt = forward_solve(model, obs.true_u);
    figure
    subplot(2,2,1)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    u = reshape(obs.true_u, [], prior.num_field);
    subplot(2,2,2)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')
    
    
    
    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = randn(prior.dof, 1);
    vmap = get_map_2020a(mlpost, mvhess, vini);
    umap = matvec_prior_L(prior, vmap) + prior.mean_u;
    solm = forward_solve(model, umap);
    
    subplot(2,2,3)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solm.state, 'edgecolor', 'none')
    u = reshape(umap, [], prior.num_field);
    subplot(2,2,4)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')
    
    
    figure
    plot(obs.data, 'o')
    hold on
    plot(solm.d, 'x')
    plot(solt.d, '-')
    
end


