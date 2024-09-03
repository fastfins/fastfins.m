%%
load_dir;

% model
h = 1/32;
a = 5; % x direction [0,5]
b = 3; % y direction [0,3]
sc = a/2; % light source center
ssigma = a/10; % light source spread
mu_s = 10;
% prior
pr_mean = -5;
pr_delta = 0.3;
pr_gamma = 1.1;

%%
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

%%
mean_u = prior.mean_u;
mesh_dof = model.mesh.dof;
delta = inv_opts.gamma;
gamma = inv_opts.cond;
sigma = inv_opts.sigma;
%
c1 = [1.5, 1.5];
c2 = [3.5, 1.5];
r1 = 1;
r2 = 0.6;
ind1 = sum((model.mesh.nodes-c1).^2, 2) < r1^2;
ind2 = sum((model.mesh.nodes-c2).^2, 2) < r2^2;
gamma = inv_opts.cond*ones(mesh_dof,1);
gamma(ind1) = 30;
gamma(ind2) = 0.2;
%
gamma = gamma.*[1, -0.99, 1];
%
prior = make_prior_mrf (model, delta, gamma, sigma);
%
prior.mean_u = mean_u;
prior.mesh_dof = mesh_dof;

%%
u = prior_random(prior, 3);
for j = 1:3
    subplot(3,1,j)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,j), 'edgecolor', 'none')
    view(2)
end

