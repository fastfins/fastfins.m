
% testing the prior, only feasible for small dof

r = randn(prior.dof, 10);
u = prior_cov_l(prior, r);
v = prior_cov_il(prior, u);
norm(v(:)-r(:))

u = prior_cov_lt(prior, r);
v = prior_cov_ilt(prior, u);
norm(v(:)-r(:))

norm(prior.L*prior.L' - prior.C)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Jd = zeros(obs.n_data, prior.dof);
for i = 1:prior.dof
    du = zeros(prior.dof, 1);
    du(i) = 1;
    Jd(:,i) = matvec_Ju(model, sol, du);
end

[U, s, V] = svd_rand_WJ(model, obs, prior, sol, 1E-10, obs.n_data);

norm( matvec_prior_Lt(prior, Jd'/obs.std)  - V*diag(s)*U')


