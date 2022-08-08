
test_case = 'tensor';
sq_flag = true;

init_func = @(x1,x2) (cos(x1.*pi.*2)-1).*sin(x2*pi*2)*2;

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

inv_opts1 = inverse_options('s2n', 10, 'cov_type', 'Conv', 'mean', 0, ...
    'centers', centers, 'radii', radii, 'weights', weights, ...
    'kernel_func', k_func, 'test_type', 'Prior');

inv_opts2 = inverse_options('s2n', 10, 'cov_type', 'MRF', 'mean', 0, ...
        'gamma', 20, 'cond', [0.2, 1, -1], 'sigma', 5);

[model, obs, prior] = setup_parabolic(model_opts, inv_opts1, test_case, []);

sol = forward_solve(model, obs.true_u);
u = reshape(obs.true_u, [], prior.num_field);
for i = 1:prior.num_field
    figure
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u(:,i), 'edgecolor', 'none')
end

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

%{
v = matvec_prior_invL(prior, obs.true_u - prior.mean_u);
[mlpt, mllkd, gmlpt, gmllkd, sol] = minus_log_post(model, obs, prior, v);

mllkd_p = zeros(prior.dof, 1);
mllkd_m = zeros(prior.dof, 1);

fd_tol = 1E-5;
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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Jdp = zeros(obs.n_data, prior.dof);
Jdm = zeros(obs.n_data, prior.dof);
tic;
for i = 1:prior.dof
    up = obs.true_u;
    up(i) = up(i)+fd_tol;
    um = obs.true_u;
    um(i) = um(i)-fd_tol;
    %
    solp = forward_solve(model, up);
    solm = forward_solve(model, um);
    Jdp(:,i) = solp.d(:);
    Jdm(:,i) = solm.d(:);
end
toc
Jd1 = (Jdp - Jdm)/(2*fd_tol);

Jd2 = zeros(obs.n_data, prior.dof);
tic;
for i = 1:prior.dof
    du = zeros(prior.dof, 1);
    du(i) = 1;
    Jd2(:,i) = matvec_Ju(model, sol, du);
end
toc

Jd3 = zeros(obs.n_data, prior.dof);
tic;
for i = 1:obs.n_data
    dy = zeros(obs.n_data, 1);
    dy(i) = 1;
    Jd3(i,:) = matvec_Jty(model, sol, dy);
end
toc

norm( Jd1(:)  - Jd2(:))
norm( Jd1(:)  - Jd3(:))
norm( Jd3(:)  - Jd2(:))

tic;
[U, s, V] = svd_rand_WJ(model, obs, prior, sol, 1E-10, obs.n_data);
toc

J1 = matvec_prior_Lt(prior,Jd1')./obs.std(:)';
J2 = matvec_prior_Lt(prior,Jd2')./obs.std(:)';
J3 = V*diag(s)*U';
norm( J1(:)  - J2(:))
%}

%{
v = matvec_prior_invL(prior, obs.true_u - prior.mean_u);
Jdp = zeros(obs.n_data, prior.dof);
Jdm = zeros(obs.n_data, prior.dof);
tic;
for i = 1:prior.dof
    up = v;
    up(i) = up(i)+fd_tol;
    um = v;
    um(i) = um(i)-fd_tol;
    %
    solp = forward_solve(model, matvec_prior_L(prior, up)+prior.mean_u);
    solm = forward_solve(model, matvec_prior_L(prior, um)+prior.mean_u);
    Jdp(:,i) = solp.d(:);
    Jdm(:,i) = solm.d(:);
end
toc
Jd1 = (Jdp - Jdm)/(2*fd_tol);

Jd2 = zeros(obs.n_data, prior.dof);
tic;
for i = 1:prior.dof
    du = zeros(prior.dof, 1);
    du(i) = 1;
    Jd2(:,i) = matvec_Ju(model, sol, matvec_prior_L(prior, du));
end
toc

tic;
[U, s, V] = svd_rand_WJ(model, obs, prior, sol, 1E-10, obs.n_data);
toc

J1 = Jd1'./obs.std(:)';
J2 = Jd2'./obs.std(:)';
J3 = V*diag(s)*U';
norm( J1(:)  - J2(:))
norm( J3(:)  - J2(:))
norm( J3(:)  - J1(:))
norm( J3(:)  - J2(:))
%}