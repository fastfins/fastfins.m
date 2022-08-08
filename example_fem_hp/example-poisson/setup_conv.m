load_dir;

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    
    test_case = 'tensor';
    sq_flag = true;
    
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
    gx  = linspace(1/4, 3/4, 4);
    gy  = linspace(1/4, 3/4, 4);
    [xx0,yy0] = meshgrid(gx, gy);
    xx0 = xx0(:);
    yy0 = yy0(:);
    
    gx  = linspace(0.75, 0.95, 2);
    gy  = linspace(0.05, 0.25, 2);
    [xx1,yy1] = meshgrid(gx, gy);
    xx1 = xx1(:);
    yy1 = yy1(:);
    
    obs_locs  = unique([xx0, yy0; xx1 yy1], 'rows');
    
    %
    model_opts = hp_options('h', 1/20, 'poly_order', 2, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);
    
    % setup the process convolution prior
    centers_x = 0.1:0.2:0.9;
    centers_y = 0.1:0.2:0.9;
    [xx,yy] = meshgrid(centers_x, centers_y);
    centers = [xx(:), yy(:)];
    radii   = ones(size(centers, 1), 1)*0.1;
    weights = ones(size(centers, 1), 1)*1.5;
    k_func  = @(x) exp( - 0.5*x.^2);
    
    inv_opts = inverse_options('s2n', 10, 'cov_type', 'Conv', 'mean', 0, ...
        'centers', centers, 'radii', radii, 'weights', weights, ...
        'kernel_func', k_func, 'test_type', 'Prior');
    
    [model, obs, prior] = setup_poisson(model_opts, inv_opts, test_case, []);
    
    solt = forward_solve(model, obs.true_u);
    figure
    subplot(2,2,1)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    u = reshape(obs.true_u, [], prior.num_field);
    for i = 1:prior.num_field
        subplot(2,2,1+i)
        trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
    end
    
    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = randn(prior.dof, 1);
    vmap = get_map_2020a(mlpost, mvhess, vini);
    umap = matvec_prior_L(prior, vmap) + prior.mean_u;
    solm = forward_solve(model, umap);
    
    figure
    subplot(2,2,1)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solm.state, 'edgecolor', 'none')
    u = reshape(umap, [], prior.num_field);
    for i = 1:prior.num_field
        subplot(2,2,1+i)
        trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
    end
    
    
    figure
    plot(obs.data, 'o')
    hold on
    plot(solm.d, 'x')
    plot(solt.d, '-')
    
end

