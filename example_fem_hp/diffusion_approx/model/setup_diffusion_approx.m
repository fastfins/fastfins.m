function [model, obs, prior] = setup_diffusion_approx(model_opts, inv_opts, output)

% Tiangang Cui, 10/Sep/2020

[p,t]   = simplemesh2d(model_opts.h, model_opts.a, model_opts.b);
model   = setup_2nd_order(p, t, model_opts);

%
model.gmres_flag  = model_opts.gmres;
model.res_tol = model_opts.res_tol;
model.explicit_ja = false;

% assemble boundary condition
model.bnd_b = zeros(model.mesh.dof, 1);
model.M_bnd = spalloc(model.mesh.dof,model.mesh.dof,0);
for bi = 1:length(model_opts.bnd_funcs)
    b_ind   = model_opts.bnd_funcs{bi}(model.mesh.nodes);
    b_force = model_opts.bc_funcs{bi}(model.mesh.nodes(:,1),model.mesh.nodes(:,2));
    if length(b_force) == 1
        b_force = ones(size(b_ind(:)))*b_force;
    else
        b_force = b_force(b_ind);
    end
    bf_weak = model.mass_bnds{bi}*sparse(b_ind(:),ones(size(b_ind(:))),b_force(:),model.mesh.dof, 1);
    model.bnd_b = model.bnd_b + bf_weak;
    model.M_bnd = model.M_bnd + model.mass_bnds{bi};
end
model.b = 2*model.bnd_b;

% apply observation operator
model.obs_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.obs_locs, model.h, false);
model.n_sensors     = size(model_opts.obs_locs,1);
model.n_datasets    = 1;

model.qoi_flag = false;

model.dof = model.mesh.dof; % for consistency
% transformation
model.exp_param = model_opts.exp_param;
model.exp_thres = model_opts.exp_thres;

model.mu_s = model_opts.mu_s;
%%%%%%%%%%%%%%%%%%% build prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prior
switch  inv_opts.cov_type
    case {'GP'}
        prior = make_prior_gp  (model.mesh, inv_opts.scale, inv_opts.power, inv_opts.sigma);
    case {'MRF'}
        prior = make_prior_mrf (model, inv_opts.gamma, inv_opts.cond, inv_opts.sigma);
    case {'Conv'}
        prior = make_prior_conv(model.mesh, inv_opts.centers, inv_opts.radii, inv_opts.weights, inv_opts.kernel_func);
end
prior.mesh_dof = model.mesh.dof;
%
if length(inv_opts.mean) == 1
    prior.mean_u = inv_opts.mean*ones(prior.mesh_dof,1);
else
    prior.mean_u = inv_opts.mean;
end

%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs.jacobian = inv_opts.jacobian;
obs.like = inv_opts.like;
%
if  ~isempty(output) && exist(output,'file')
    tmp         = load(output);
    obs.data    = tmp.data;
    obs.std     = tmp.std;
    obs.true_u  = tmp.true_u;
else
    r = randn(prior.dof, 1);
    true_u  = matvec_prior_L(prior, r) + prior.mean_u;
    soln    = forward_solve(model, true_u);  % reference solution
    
    % the s.t.d. is calculated from the signal to noise ratio
    if inv_opts.s2n > 0
        std = mean(abs(soln.d(:)))/inv_opts.s2n;
    else
        std = inv_opts.std;
    end
    
    % generate data
    data    = soln.d + randn(model.n_sensors, model.n_datasets)*std;
    
    obs.data    = data;
    obs.std     = std;
    obs.true_u  = true_u;
end

if  ~isempty(output) && ~exist(output,'file')
    save(output, 'data', 'std', 'true_u');
end

obs.n_sensors    = model.n_sensors;
obs.n_datasets   = model.n_datasets;
obs.n_data       = obs.n_sensors*obs.n_datasets;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

