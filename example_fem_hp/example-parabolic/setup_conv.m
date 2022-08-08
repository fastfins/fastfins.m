load_dir;

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    
    test_case = 'tensor';
    sq_flag = true;
    
    init_func = @(x1,x2) (cos(x1.*pi.*2)-1).*sin(x2*pi*2)*2;
    
    % boundary conditions
    %{
    tol = 1E-4;
    bnd_funcs = {@(p) find(abs(p(:,2))<tol), @(p) find(abs(p(:,2)-1)<tol)};
    bc_types  = {'essential', 'essential'};
    bc_funcs  = {@(x1,x2)-(cos(x1.*pi.*2)-1)/sqrt(2), @(x1,x2)cos(x1.*pi.*2)-1};
    %}
    
    %
    %{
    f1  = @(x1,x2) exp(-0.5*((x1-0.8).^2+(x2-0.5).^2)/0.05^2)'/(2*pi*0.05^2);
    f2  = @(x1,x2) exp(-0.5*((x1-0.2).^2+(x2-0.5).^2)/0.05^2)'/(2*pi*0.05^2);
    for_func  = @(x1,x2) f1(x1,x2) - f2(x1,x2);
    %}
    
    % flux function for qoi
    flux_fx = @(x) uniform_fun (x, 0, 1);
    flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
    flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);
    
    % obs locations
    gx  = linspace(1/4, 3/4, 2);
    gy  = linspace(1/4, 3/4, 2);
    [xx0,yy0] = meshgrid(gx, gy);
    xx0 = xx0(:);
    yy0 = yy0(:);
       
    obs_locs  = unique([xx0, yy0; 0.8, 0.2], 'rows');
    
    dh = 1/20;
    dt = 1/100;
    %
    %{
    model_opts = hp_options('h', dh, 'poly_order', 2, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs, ...
        'time_dep', true, 'dt', dt, 'init_func', init_func, 'force_func', for_func, ...
        'obs_tstart', 0.05, 'obs_tfinal', 0.15, 'obs_ntime', 6, 'pred_t', 0.2);
    %}
    model_opts = hp_options('h', dh, 'poly_order', 2, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'time_dep', true, 'dt', dt, 'init_func', init_func, ...
        'obs_tstart', 0.05, 'obs_tfinal', 0.15, 'obs_ntime', 6, 'pred_t', 0.2);
    
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
    
    [model, obs, prior] = setup_parabolic(model_opts, inv_opts, test_case, []);
    
    solt = forward_solve(model, obs.true_u);
    figure(1)
    u = reshape(obs.true_u, [], prior.num_field);
    for i = 1:prior.num_field
        subplot(2,2,i)
        trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
        view(3)
    end
    figure(2)
    for i = 1:model.T_nstep
        figure(2)
        trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state(:,i), 'edgecolor', 'none')
        view(2)
        title(['t=' num2str(model.dt*i)])
        pause(0.1)
    end
    
    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = randn(prior.dof, 1);
    vmap = get_map_2020a(mlpost, mvhess, vini);
    umap = matvec_prior_L(prior, vmap) + prior.mean_u;
    solm = forward_solve(model, umap);
    
    figure(3)
    u = reshape(umap, [], prior.num_field);
    for i = 1:prior.num_field
        subplot(2,2,i)
        trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
        view(3)
    end
    figure(4)
    for i = 1:model.T_nstep
        figure(4)
        trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solm.state(:,i), 'edgecolor', 'none')
        view(2)
        title(['t=' num2str(model.dt*i)])
        pause(0.1)
    end

    
    figure
    plot(obs.data, 'o')
    hold on
    plot(solm.d, 'x')
    plot(solt.d, '-')
    
end

