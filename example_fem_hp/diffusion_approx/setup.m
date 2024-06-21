%%
load_dir;

% model
hs = 1./[4, 8, 16, 32, 64];
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
for i = 1:length(hs)
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
    model_opts = hp_options('h', hs(i), 'poly_order', 1, 'quad_order', 15, ...
        'xyratio', 1, 'obs_locs', obs_locs, ...
        'bnd_funcs', bnd_funcs, 'bc_funcs', bc_funcs);

    inv_opts = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', -5, ...
        'gamma', pr_delta, 'cond', pr_gamma, 'sigma', 1);

    model_opts.mu_s = mu_s;
    model_opts.a = a;
    model_opts.b = b;

    [models{i}, obss{i}, priors{i}] = setup_diffusion_approx(model_opts, inv_opts, []);

end

%%
for i = 1:length(hs)
    solt = forward_solve(models{i}, obss{i}.true_u);
    figure
    subplot(1,2,1)
    trisurf(models{i}.mesh.node_tri, models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    u = reshape(obss{i}.true_u, [], 1);
    subplot(1,2,2)
    trisurf(models{i}.mesh.node_tri, models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2), u, 'edgecolor', 'none')
end

%%
% access the mesh of the first model
i = 1;
models{i}.mesh % the mesh structure
models{i}.mesh.nodes; % num_nodes x 2, coordinates of nodes
models{i}.mesh.node_map; % num_elem x 3, node indices for each element
j = 1;
% x and y coordindates of the nodes of jth element
xj0 = models{i}.mesh.nodes(models{i}.mesh.node_map(j,:), 1);
yj0 = models{i}.mesh.nodes(models{i}.mesh.node_map(j,:), 2);
% x and y coordindates of the nodes of j+1th element
xj1 = models{i}.mesh.nodes(models{i}.mesh.node_map(j+1,:), 1);
yj1 = models{i}.mesh.nodes(models{i}.mesh.node_map(j+1,:), 2);

figure
plot(xj0, yj0, 'ko-', 'LineWidth', 3)
hold on
plot(xj1, yj1, 'bo:', 'LineWidth', 3)
title('counter clock node ordering, also showing the cut of the square')
