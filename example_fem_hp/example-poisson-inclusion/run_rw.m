load obs

ml_target = @(v) minus_log_post(model, obs, prior_inclusion, v);
%
stat.n = 1E2;
stat.cross = eye(prior_inclusion.dof)*stat.n;
stat.sum = u_true(:)*stat.n;
%
opt = mcmc_options('proposal', 'RW', 'nstep', 1E6, 'sigma', -4);

tic;
[out, stat] = amcmc(ml_target, stat, u_true(:), opt);
toc

ind = 1E5:5E2:1E6;
kk = inclusions(model, out.samples(:,ind));
k_mean = mean(kk, 2);
k_std = std(kk, [], 2);


figure
subplot(2,2,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), k_true, 'edgecolor', 'none')
view(2)

subplot(2,2,3)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), k_mean, 'edgecolor', 'none')
view(2)

subplot(2,2,4)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), k_std, 'edgecolor', 'none')
view(2)

