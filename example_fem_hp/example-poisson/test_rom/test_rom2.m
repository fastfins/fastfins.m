

test_case = 'scalar';
sq_flag = false;
scale_used = 30; 
h_opt = 1/32; 

% boundary conditions
tol = 1E-5;

bnd_funcs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,1)-1)<tol)};
bc_types  = {'essential', 'essential'};
%bc_funcs = {@(x1,x2) 1+x2/2, @(x1,x2) -sin(2*pi*x2)-1};
bc_funcs  = {@(x1,x2)-(cos(x2.*pi.*2 - pi/2)-1)/sqrt(2), @(x1,x2)cos(x2.*pi.*2 + pi/2)-1};
%bc_funcs = {@(x1,x2)-((cos(x2.*pi.*4)-1).*(x2>0.5))/sqrt(2), @(x1,x2)(cos(x2.*pi.*4)-1).*(x2<0.5)};

% flux function for qoi
flux_fx = @(x) uniform_fun (x, 0, 1);
flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);

% obs locations
%{
% old one
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
%}

% single point
obs_locs = [0.5,0.5]; 

%
model_opts = hp_options('h', h_opt, 'poly_order', 2, 'quad_order', 15, ...
    'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
    'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);


% TC's old test prior
%{
% setup the process convolution prior
centers_x = 0.2:0.3:0.8;
centers_y = 0.2:0.3:0.8;
[xx,yy] = meshgrid(centers_x, centers_y);
centers = [xx(:), yy(:)];
radii   = ones(size(centers, 1), 1)*0.3;
if sq_flag
    weights = ones(size(centers, 1), 1)*0.5;
else
    weights = ones(size(centers, 1), 1)*1;
end
k_func  = @(x) exp( - 0.5*x.^2);
inv_opts = inverse_options('s2n', 10, 'cov_type', 'Conv', 'mean', 0, ...
    'centers', centers, 'radii', radii, 'weights', weights, ...
    'kernel_func', k_func, 'test_type', 'Prior');
%}

% new prior
inv_opts = inverse_options('s2n', 10, 'cov_type', 'GP', 'mean', 0, ...
    'scale', scale_used, 'power', 2, 'sigma', 1);
%
[model, true_obs, prior] = setup_poisson(model_opts, inv_opts, test_case, []);
prior_KL = basis_KL(prior, 1-1E-2);

solt = forward_solve(model, true_obs.true_u);
figure
subplot(2,2,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
u = reshape(true_obs.true_u, [], prior.num_field);
for i = 1:prior.num_field
    subplot(2,2,1+i)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
end

% run the map using prior_KL
mlpost  = @(v) minus_log_post(model, true_obs, prior_KL, v);
mvhess  = @(sol, dv) matvec_PPFisher(model, prior_KL, sol, dv);
vini = randn(prior_KL.dof, 1);
vmap = get_map_2020a(mlpost, mvhess, vini);
umap = matvec_prior_L(prior_KL, vmap) + prior_KL.mean_u;
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
plot(true_obs.data, 'o')
hold on
plot(solm.d, 'x')
plot(solt.d, '-')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1E3;
vs = randn(prior_KL.dof, n);
%
%
states = zeros(model.dof,n);
disp(['evaluate ' num2str(n) ' full models for building rom'])
tic;
for i = 1:n
    u   = matvec_prior_L(prior_KL, vs(:,i)) + prior_KL.mean_u;
    % forward solve
    sol = forward_solve(model, u);
    states(:,i) = sol.state;
end
toc
weights = ones(1,n)/n;
%
prior_redu.bases{1} = prior_KL.basis;
prior_redu.mean_us{1} = prior_KL.mean_u;
%
sub_vs{1} = vs;
%

% reduce these to increase accuracy
rom_opts.deim_reg_factor = 1.3;
rom_opts.prior_deim_tol  = 1E-6;
rom_opts.pod_state_tol   = 1E-4;
disp('build rom')
tic;
rom = setup_poisson_rom(model, prior_redu, sub_vs, states, weights, test_case, rom_opts);
toc
disp(['rom state dim = ' num2str(rom.dof) ', rom DEIM dim = ' num2str(size(rom.Ks{1},1))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_d = zeros(true_obs.n_data, n);
rom_d  = zeros(true_obs.n_data, n);
vs_new = randn(prior_KL.dof, n);
disp('run full models')
tic;
for i = 1:n
    u   = matvec_prior_L(prior_KL, vs_new(:,i)) + prior_KL.mean_u;
    % forward solve
    sol = forward_solve(model, u);
    full_d(:,i) = sol.d;
end
toc
disp('run reduced models')
tic;
for i = 1:n
    solr = rom_solve(rom, vs_new(:,i));
    rom_d(:,i) = solr.d;
end
toc

figure
plot(full_d(:) - rom_d(:))
hold on
plot(full_d(:))
plot(rom_d(:))
legend('full - rom on observables', 'full', 'rom')

disp(['relative error:' num2str(norm(full_d(:) - rom_d(:))/norm(rom_d(:)))])
