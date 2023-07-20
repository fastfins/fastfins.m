load_dir;

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    
    test_case = 'scalar';
    sq_flag = false;

    % boundary conditions
    tol = 1E-5;
    %bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
    %bc_types  = {'essential', 'essential'};
    %bc_funcs  = {@(x1,x2) 1, @(x1,x2) -1};

    bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
    bc_types  = {'essential', 'essential'};
    bc_funcs  = {@(x1,x2)-(cos(x2.*pi.*2)-1)/sqrt(2), @(x1,x2)cos(x2.*pi.*2)-1};
    
    % flux function for qoi
    flux_fx = @(x) uniform_fun (x, 0, 1);
    flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
    flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);
    
    % obs locations
    gx  = linspace(0, 1, 9);
    gy  = linspace(0, 1, 9);
    [xx0,yy0] = meshgrid(gx(2:end-1), gy(2:end-1));
    xx0 = xx0(:);
    yy0 = yy0(:);
    
    obs_locs  = unique([xx0, yy0], 'rows');
    
    %
    model_opts = hp_options('h', 1/16, 'poly_order', 2, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);
    
    %inv_opts = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', 0, ...
    %    'gamma', 20, 'cond', [1, 1, 0], 'sigma', 1);

    inv_opts = inverse_options('s2n', 10, 'cov_type', 'GP', 'mean', 0, ...
        'scale', 50, 'power', 2, 'sigma', 2);
    
    [model, obs, prior] = setup_poisson(model_opts, inv_opts, test_case, []);
    
    solt = forward_solve(model, obs.true_u);
    
    figure
    subplot(2,2,1)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    subplot(2,2,2)
    u = reshape(obs.true_u, [], prior.num_field);
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')


    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = randn(prior.dof, 1);
    vmap = get_map_2023a(mlpost, mvhess, vini);
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

