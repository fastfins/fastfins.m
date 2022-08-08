
%
fin_path = /Users/tcui001/Work/fastfins.m';
warning('off')
rmpath(genpath(fin_path));
warning('on')
addpath(genpath([fin_path '/fastfins/']));
addpath(genpath([fin_path '/fem_hp']));
%
addpath(genpath([fin_path '/example_fem_hp/example-poisson']));
addpath(genpath([fin_path '/rom_tools']));

fastfins_check_solvers()

renew_setup = true;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    
    test_case = 'scalar'; % 'tensor'; % 
    sq_flag = false;
    
    % boundary conditions
    tol = 1E-5;
    bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
    bc_types  = {'essential', 'essential'};
    bc_funcs  = {@(x1,x2)-(cos(x2.*pi.*2)-1)/sqrt(2), @(x1,x2)cos(x2.*pi.*2)-1};
    

    %bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
    %bc_types  = {'essential', 'essential'};
    %bc_funcs  = {@(x1,x2) 1, @(x1,x2) 0};

    % flux function for qoi
    flux_fx = @(x) uniform_fun (x, 0, 1);
    flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
    flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);
    
    % obs locations
    gx  = linspace(0, 1, 17);
    gy  = linspace(0, 1, 17);
    [xx0,yy0] = meshgrid(gx(2:end-1), gy(2:end-1));
    xx0 = xx0(:);
    yy0 = yy0(:);
    
    obs_locs  = unique([xx0, yy0], 'rows');
    
    %
    model_opts = hp_options('h', 1/32, 'poly_order', 2, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);
    
    % setup the process convolution prior
    ngrid = 10;
    centers_x = (1/ngrid/2):(1/ngrid):(1-1/ngrid/2);
    centers_y = (1/ngrid/2):(1/ngrid):(1-1/ngrid/2);
    [xx,yy] = meshgrid(centers_x, centers_y);
    centers = [xx(:), yy(:)];
    radii   = ones(size(centers, 1), 1)*0.05;
    weights = ones(size(centers, 1), 1)*2;
    k_func  = @(x) exp( - 0.5*x.^2);
    
    inv_opts = inverse_options('s2n', 20, 'cov_type', 'Conv', 'mean', 0, ...
        'centers', centers, 'radii', radii, 'weights', weights, ...
        'kernel_func', k_func, 'test_type', 'Prior');
    
    [model, obs, prior] = setup_poisson(model_opts, inv_opts, test_case, []);
    % drill a channel
    channel_u = image_channel(model.mesh, 1E-2, 1E2, 0.075);
    channel_u = log(channel_u);
    % new data
    solt = forward_solve(model, channel_u);  % reference solution
    % the s.t.d. is calculated from the signal to noise ratio
    std = mean(abs(solt.d(:)))/inv_opts.s2n; %(1-xx0)/inv_opts.s2n;
    % generate data
    data = solt.d + randn(model.n_sensors, model.n_datasets).*std;
    %
    obs.data    = data;
    obs.std     = std;
    obs.true_u  = channel_u;
    % end of channel setup

    figure
    subplot(2,2,1)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    u = reshape(obs.true_u, [], prior.num_field);
    subplot(2,2,2)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')

    mlpost  = @(v) minus_log_post(model, obs, prior, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior, sol, dv);
    vini = randn(prior.dof, 1); % prior.basis\obs.true_u;% 
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

