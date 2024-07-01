%%
load_dir;

% model
h = 1./8;
a = 5; % x direction [0,5]
b = 3; % y direction [0,3]
sc = a/2; % light source center
ssigma = a/10; % light source spread
mu_s = 10;
% prior
pr_mean = -5;
pr_delta = 0.3;
pr_gamma = 1.1;

% boundary conditions
tol = 1E-4;
bnd_funcs = {
    @(p) find(abs(p(:,2))<tol), ...     % y = 0, bottom
    @(p) find(abs(p(:,1)-a)<tol), ...   % x = a, right
    @(p) find(abs(p(:,2)-b)<tol), ...   % y = b, top
    @(p) find(abs(p(:,1))<tol)};        % x = 0, left
bc_funcs  = {
    @(x1,x2) 0, ... % y = 0, bottom
    @(x1,x2) 0, ... % x = a, right
    @(x1,x2) 6*exp(-0.5*(x1-sc).^2/ssigma^2), ... % y = b, top
    @(x1,x2) 0};    % x = 0, left

% obs locations, these are only for debugging purpose, will be replaced
% later on
gx  = linspace(1/4, 3/4, 4)*a;
gy  = linspace(1/4, 3/4, 4)*b;
[xx0,yy0] = meshgrid(gx, gy);
xx0 = xx0(:);
yy0 = yy0(:);

gx  = linspace(0.75, 0.95, 2)*a;
gy  = linspace(0.05, 0.25, 2)*b;
[xx1,yy1] = meshgrid(gx, gy);
xx1 = xx1(:);
yy1 = yy1(:);

obs_locs  = unique([xx0, yy0; xx1 yy1], 'rows');

% piecewise linear
model_opts = hp_options('h', h, 'poly_order', 1, 'quad_order', 15, ...
    'xyratio', 1, 'obs_locs', obs_locs, ...
    'bnd_funcs', bnd_funcs, 'bc_funcs', bc_funcs);

inv_opts = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', -5, ...
    'gamma', pr_delta, 'cond', pr_gamma, 'sigma', 1);

model_opts.mu_s = mu_s;
model_opts.a = a;
model_opts.b = b;

[model, obs, prior] = setup_diffusion_approx(model_opts, inv_opts, []);

solt = forward_solve(model, obs.true_u);
figure
subplot(1,2,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
u = reshape(obs.true_u, [], 1);
subplot(1,2,2)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')

%%
prior_KL = basis_KL(prior, 1-1E-2);
N = 1E3;
sub_vs = randn(prior_KL.dof, N);
weights = ones(N, 1);
states = zeros(model.dof, N);
ds = zeros(obs.n_data, N);
tic
for i = 1:N
    sol = forward_solve(model, prior_KL.basis*sub_vs(:,i)+prior_KL.mean_u);
    states(:,i) = sol.state;
    ds(:,i) = sol.d;
end
toc

%%
rom_opts.pod_state_tol = 1E-4;
rom_opts.prior_deim_tol = 1E-6;
rom_opts.deim_reg_factor = 2;
rom = setup_rom(model, states, prior_KL, sub_vs, weights, rom_opts);

%%
rstates = zeros(rom.dof, N);
rds = zeros(obs.n_data, N);
tic
for i = 1:N
    rsol = rom_solve(rom, sub_vs(:,i));
    %rsol = rom_solve_debug(model, prior_KL, rom, sub_vs(:,i));
    rstates(:,i) = rsol.state;
    rds(:,i) = rsol.d;
end
toc
