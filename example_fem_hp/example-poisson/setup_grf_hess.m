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
    gx  = linspace(0, 5/8, 11);
    gy  = linspace(0, 5/8, 11);
    [xx0,yy0] = meshgrid(gx(2:end-1), gy(2:end-1));
    xx0 = xx0(:);
    yy0 = yy0(:);
    
    obs_locs  = unique([xx0, yy0], 'rows');
    
    %
    model_opts = hp_options('h', 1/16, 'poly_order', 1, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);
    
    %inv_opts = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', 0, ...
    %    'gamma', 20, 'cond', [1, 1, 0], 'sigma', 1);

    inv_opts = inverse_options('s2n', 10, 'cov_type', 'GP', 'mean', 0, ...
        'scale', 50, 'power', 1, 'sigma', 2);
    
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

%%
tol = 1E-5;
Jd1 = zeros(numel(vmap)*obs.n_data,numel(vmap));
Jd2 = zeros(numel(vmap)*obs.n_data,numel(vmap));
for i = 1:numel(vmap)
    vp = vmap;
    vp(i) = vp(i) + tol;
    vm = vmap;
    vm(i) = vm(i) - tol;
    [~,~,~,~,solp] = minus_log_post(model, obs, prior, vp);
    Jp = explicit_jacobian(model, solp);
    [~,~,~,~,solm] = minus_log_post(model, obs, prior, vm);
    Jm = explicit_jacobian(model, solm);
    %
    tmp = matvec_prior_Lt(prior, (Jp-Jm)'/(2*tol))';
    ind = (1:obs.n_data) + obs.n_data*(i-1);
    Jd1(ind,:) = tmp;
    Jd2(:,i) = reshape(tmp, [], 1);
    %Jd(ind,:) = (Jp-Jm)/(2*tol);
end

%%
es = eye(numel(vv),10);


%%
r = randn(size(vmap,1),10);
u = matvec_prior_L(prior, r);
[~,~,~,~,sol] = minus_log_post(model, obs, prior, vmap);
%Ju = explicit_jacobian(model, sol);
%Ju = Ju./obs.std(:);
%Jv = matvec_prior_Lt(prior, Ju')';
%
a = matvec_PPJ2(model,prior,sol,r);
for i = 1:size(r,2)
    ad = reshape(Jd2*r(:,i), obs.n_data, []);

    t0 = ad'*sol.gd;
    t1 = a{i}'*sol.gd;
    t2 = matvec_dyJ2du(model, sol, sol.gd, u(:,i));
    t2 = matvec_prior_Lt(prior, t2); 
    [norm(t0-t1)/norm(t0) norm(t1-t2)/norm(t0)]
end

%%
vv = vmap + randn(size(vmap))*1E-2;
H = zeros(numel(vv));
for i = 1:numel(vv)
    vp = vv;
    vp(i) = vp(i) + tol;
    vm = vv;
    vm(i) = vm(i) - tol;
    [~,~,~,gmllkd_p] = minus_log_post(model, obs, prior, vp);
    [~,~,~,gmllkd_m] = minus_log_post(model, obs, prior, vm);
    H(:,i) = (gmllkd_p-gmllkd_m)/(2*tol);
end
H = (H+H')/2;

%%
HH = zeros(numel(vv));
[~,~,~,~,sol] = minus_log_post(model, obs, prior, vv);
for i = 1:numel(vv)
    e = zeros(size(vv));
    e(i) = 1;
    HH(:,i) = matvec_PPHess(model, prior, sol, e);
end
HH = (HH+HH')/2;

norm(HH-H)

%%
nmax = 50;
%
[V1,d1] = eigen_PPHess_FD(model, obs, prior, vmap, 1E-6, 1E-10, nmax);
[V2,d2] = eigen_PPHess   (model, obs, prior, vmap, 1E-10, nmax);
%
figure
plot(d1, 'x'), hold on, plot(d2, 'o')
figure
plot((V1(:,1:20)'*V2(:,1:20)))