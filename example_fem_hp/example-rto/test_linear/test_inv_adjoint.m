
sq_flag = true;

% test case
test_case = 'tensor';


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
model_opts = hp_options('h', 1/10, 'poly_order', 2, 'quad_order', 15, ...
    'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
    'bnd_funcs', bnd_funcs,'bc_types', bc_types, 'bc_funcs', bc_funcs);

inv_opts = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', 0, ...
        'gamma', 20, 'cond', [1, 0.5, -0.5], 'sigma', 2);

[model, obs, prior] = setup_poisson(model_opts, inv_opts, test_case, []);

sol = forward_solve(model, obs.true_u);
u = reshape(obs.true_u, [], prior.num_field);
for i = 1:prior.num_field
figure
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
end
figure
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), sol.state, 'edgecolor', 'none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% testing the prior, only feasible for small dof

r = randn(prior.dof, 10);
u = prior_cov_l(prior, r);
v = prior_cov_il(prior, u);
norm(v(:)-r(:))

u = prior_cov_lt(prior, r);
v = prior_cov_ilt(prior, u);
norm(v(:)-r(:))

L = prior_cov_l(prior, eye(prior.dof));
Lt = prior_cov_lt(prior, eye(prior.dof));
Q = prior.K*(prior.M\prior.K);

norm(L'-Lt)
I = Lt*Q*L; % this needs to be near identity % pass the test
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = matvec_prior_invL(prior, obs.true_u - prior.mean_u);
[mlpt, mllkd, gmlpt, gmllkd, sol] = minus_log_post(model, obs, prior, v);

mllkd_p = zeros(prior.dof, 1);
mllkd_m = zeros(prior.dof, 1);

fd_tol = 1E-6;
tic;
for i = 1:prior.dof
    vp = v;
    vp(i) = vp(i)+fd_tol;
    vm = v;
    vm(i) = vm(i)-fd_tol;
    [~,mllkd_p(i)] = minus_log_post(model, obs, prior, vp);
    [~,mllkd_m(i)] = minus_log_post(model, obs, prior, vm);
end
toc

fd_gmllkd = (mllkd_p - mllkd_m)/(2*fd_tol);
norm(fd_gmllkd - gmllkd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


J = explicit_jacobian(model, sol);

Jd = zeros(obs.n_data, prior.dof);
tic;
for i = 1:prior.dof
    du = zeros(prior.dof, 1);
    du(i) = 1;
    Jd(:,i) = matvec_Ju(model, sol, du);
end
toc

norm(J - Jd)

[U, s, V] = svd_rand_WJ(model, obs, prior, sol, 1E-10, obs.n_data);

norm( matvec_prior_Lt(prior, J'/obs.std)  - V*diag(s)*U')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = matvec_prior_invL(prior, obs.true_u - prior.mean_u);
[~, ~, ~, ~, sol] = minus_log_post(model, obs, prior, v);

c = sol.qoi*2;
sigma = 0.5;
[~, gmlr] = minus_log_rare(model, prior, c, sigma, v);
[~, ~, ~, gmllkd] = minus_log_post_rare(model, obs, prior, c, sigma, v);

mlr_p = zeros(prior.dof, 1);
mlr_m = zeros(prior.dof, 1);
mllkd_p = zeros(prior.dof, 1);
mllkd_m = zeros(prior.dof, 1);

fd_tol = 1E-5;
tic;
for i = 1:prior.dof
    vp = v;
    vp(i) = vp(i)+fd_tol;
    vm = v;
    vm(i) = vm(i)-fd_tol;
    
    mlr_p(i) = minus_log_rare(model, prior, c, sigma, vp);
    mlr_m(i) = minus_log_rare(model, prior, c, sigma, vm);
    
    [~, mllkd_p(i)] = minus_log_post_rare(model, obs, prior, c, sigma, vp);
    [~, mllkd_m(i)] = minus_log_post_rare(model, obs, prior, c, sigma, vm);

end
toc

fd_gmlr = (mlr_p - mlr_m)/(2*fd_tol);
norm(fd_gmlr - gmlr)

fd_gmllkd = (mllkd_p - mllkd_m)/(2*fd_tol);
norm(fd_gmllkd - gmllkd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vs = v*10.^linspace(-5, 1, 50);
for i = 1:size(vs, 2)
    [mlrs(i), ~, sol] = minus_log_rare(model, prior, c, sigma, vs(:,i));
    qs(i) = sol.qoi;
end
