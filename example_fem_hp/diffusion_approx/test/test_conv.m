%%

% model
hs = 1./[4, 8, 16, 32, 64];
a = 5; % x direction [0,5]
b = 3; % y direction [0,3]
% prior
pr_mean = -5;
pr_delta = 0.3;
pr_gamma = 1.1;

C = 10;
f = @(x,y) sqrt(x+y+C);
true_sol = @(x,y) exp(-f(x,y)); 
true_mua = @(x,y) 1./(2*f(x,y));
mu_s = 0;

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
        @(x1,x2) (1/pi+1/4)*true_sol(x1,0), ... % y = 0, bottom
        @(x1,x2) (1/pi-1/4)*true_sol(a,x2), ... % x = a, right
        @(x1,x2) (1/pi-1/4)*true_sol(x1,b), ... % y = b, top
        @(x1,x2) (1/pi+1/4)*true_sol(0,x2)};    % x = 0, left

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

    [models{i},obss{i},priors{i}] = setup_diffusion_approx(model_opts, inv_opts, []);

end

%%
es1 = [];
es2 = [];
for i = 1:length(hs)
    true_u = true_mua(models{i}.mesh.nodes(:,1),models{i}.mesh.nodes(:,2));
    true_s = true_sol(models{i}.mesh.nodes(:,1),models{i}.mesh.nodes(:,2));
    solt = forward_solve(models{i}, log(true_u));
    diff = solt.state - true_s;
    figure
    subplot(2,2,1)
    trisurf(models{i}.mesh.node_tri, models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    title('numerical');
    subplot(2,2,2)
    trisurf(models{i}.mesh.node_tri, models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2), true_u, 'edgecolor', 'none')
    title('\mu_a');
    subplot(2,2,3)
    trisurf(models{i}.mesh.node_tri, models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2), true_s, 'edgecolor', 'none')
    title('true');
    subplot(2,2,4)
    trisurf(models{i}.mesh.node_tri, models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2), diff, 'edgecolor', 'none')
    title('diff');

    es1(i) = sqrt(diff'*models{i}.mass*diff);
    es2(i) = norm(models{i}.obs_operator*(true_u.*diff));
end

figure
loglog(hs, es1, 'x-')
hold on
loglog(hs, es2, 'o-')
title('L^2 error vs h')
legend('state', 'observable')
