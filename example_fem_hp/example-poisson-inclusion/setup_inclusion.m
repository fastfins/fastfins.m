setup_rf;

model.inclusion_flag = true;
model.scale = [1,1,0.5,1,1,0.5,4,1]'; % center, r, center, r, kin, kout
model.shift = [0,0,0.1,0,0,0.1,1,1E-5]';

d = 8;
prior_inclusion.dof = d;
prior_inclusion.type = 'Field';
prior_inclusion.cov_type = 'GP';
prior_inclusion.L = eye(d);
prior_inclusion.C = eye(d);
prior_inclusion.mean_u = zeros(d,1);

samples = prior_random(prior_inclusion, 9);
k = inclusions(model, samples);
    
figure
for i = 1:9
    subplot(3,3,i)
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), k(:,i), 'edgecolor', 'none')
    view(2)
end

u_true = [0.5, -0.2, -0.5, -0.6, -0.6, -1, 3, 0.1]';
k_true = inclusions(model, u_true);
figure
subplot(2,2,1)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), k_true, 'edgecolor', 'none')
view(2)

sol = forward_solve(model, u_true);  % reference solution

%{
s2n = 20;
std = mean(abs(sol.d(:)))/s2n;
% generate data
data = sol.d + randn(model.n_sensors, model.n_datasets)*std;
obs.data = data;
obs.std = std;
obs.true_u = u_true;
%}

subplot(2,2,2)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), sol.kappa, 'edgecolor', 'none')
subplot(2,2,3)
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), sol.state, 'edgecolor', 'none')

%[out, stat] = amcmc(ml_target, stat, init, options);