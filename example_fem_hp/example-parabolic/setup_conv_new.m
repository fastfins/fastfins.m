load_dir;

%%
test_case = 'scalar';
sq_flag = false;

%init_func = @(x1,x2) (cos(x1.*pi.*2)-1).*sin(x2*pi*2)*2;

sigma = 0.1;
init_func = @(x1,x2) exp( -0.5*((x1-0.3).^2+(x2-0.3).^2)/sigma^2 )*2 ...
    + exp( -0.5*((x1-0.7).^2+(x2-2.7).^2)/sigma^2 ) ...
    - exp( -0.5*((x1-0.5).^2+(x2-1.5).^2)/sigma^2 );

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
flux_fy = @(y) triangle_fun(y, 0, 3, 'right');
flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);

% obs locations
gx1 = linspace(1/5, 4/5, 4);
gy1 = linspace(1/5, 4/5, 4);
[xx1,yy1] = meshgrid(gx1, gy1);
xx1 = xx1(:);
yy1 = yy1(:);
%
gx2 = linspace(1/5, 4/5, 2);
gy2 = linspace(1/5, 4/5, 2) + 2;
[xx2,yy2] = meshgrid(gx2, gy2);
xx2 = xx2(:);
yy2 = yy2(:);

obs_locs = unique([xx1, yy1; xx2, yy2], 'rows');
pred_locs = unique([0.5, 0.5; 0.5, 1.5; 0.5, 2.5], 'rows');

dh = 1/32;
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
    'xyratio', 3, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
    'time_dep', true, 'dt', dt, 'init_func', init_func, ...
    'obs_tstart', 0.05, 'obs_tfinal', 0.15, 'obs_ntime', 6, 'pred_t', 0.25);

model_opts.pred_locs = pred_locs;

% setup the process convolution prior
centers_x = 0.1:0.1:0.9;
centers_y = 0.1:0.1:2.9;
[xx,yy] = meshgrid(centers_x, centers_y);
centers = [xx(:), yy(:)];
radii   = ones(size(centers, 1), 1)*0.1;
weights = ones(size(centers, 1), 1)*1;
k_func  = @(x) exp(-0.5*x.^2);

inv_opts = inverse_options('s2n', 10, 'cov_type', 'Conv', 'mean', 0, ...
    'centers', centers, 'radii', radii, 'weights', weights, ...
    'kernel_func', k_func, 'test_type', 'Prior');

[model, obs, prior] = setup_parabolic(model_opts, inv_opts, test_case, []);

%%
tic;
solt = forward_solve(model, obs.true_u);
toc

%%
figure(1)
n = 4;
u = prior_random(prior, n);
for i = 1:n
    subplot(n,1,i)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), u(:,i), 'edgecolor', 'none')
    view(2)
end

%%
f2 = figure('Position',[100, 100, 800, 800]);
subplot(3,1,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), obs.true_u, 'edgecolor', 'none')
subplot(3,1,2)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), model.init, 'edgecolor', 'none')
subplot(3,1,3)
for i = 1:model.T_nstep
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), solt.state(:,i), 'edgecolor', 'none')
    %view(2)
    title(['t=' num2str(model.dt*i)])
    pause(0.1)
end

%%
f3 = figure('Position',[100, 800, 800, 800]);
subplot(3,1,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), obs.true_u, 'edgecolor', 'none')
hold on
plot3(model_opts.obs_locs(:,2), model_opts.obs_locs(:,1), max(obs.true_u)*ones(size(model_opts.obs_locs,1),1), 'r.', markersize=20)
plot3(model_opts.pred_locs(:,2), model_opts.pred_locs(:,1), max(obs.true_u)*ones(size(model_opts.pred_locs,1),1), 'rx', markersize=10)
view(2)
subplot(3,1,2)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), model.init, 'edgecolor', 'none')
hold on
plot3(model_opts.obs_locs(:,2), model_opts.obs_locs(:,1), max(obs.true_u)*ones(size(model_opts.obs_locs,1),1), 'r.', markersize=20)
plot3(model_opts.pred_locs(:,2), model_opts.pred_locs(:,1), max(obs.true_u)*ones(size(model_opts.pred_locs,1),1), 'rx', markersize=10)
view(2)
subplot(3,1,3)
for i = 1:model.T_nstep
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,2), model.mesh.nodes(:,1), solt.state(:,i), 'edgecolor', 'none')
    view(2)
    %view(2)
    title(['t=' num2str(model.dt*i)])
    pause(0.1)
end

%%
figure
subplot(1,2,1)
plot(model.obs_ind*model.dt, reshape(obs.data, [], numel(model.obs_ind)), 'o')
hold on
plot(model.obs_ind*model.dt, reshape(solt.d, [], numel(model.obs_ind)), '-x')
subplot(1,2,2)
plot(model.pred_ind*model.dt, reshape(solt.q, [], numel(model.pred_ind)), '-x')
