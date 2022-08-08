function[model, obs, prior] = setup_p_poisson(model_opts, inv_opts, p_rate, epsilon, output)
%setup_model
%
% setup data structures used for model computation
%
% Tiangang Cui, 10/Sep/2020

[p,t] = simplemesh2d(model_opts.h, 1, model_opts.xyratio);
model = setup_2nd_order(p, t, model_opts);
%
model.WdetJ = model.local_elem.quad_weights(:)*model.mesh.detJ(:)';
%
model.pa_bnd_i      = 1; % must be specified (not in the options)
model.pa_bnd_mesh   = model.bnd_mesh{model.pa_bnd_i};
%
model.pa_bnd_WdetJ  = model.local_elem_bnd.quad_weights(:) * model.pa_bnd_mesh.detJ(:)';
%
model.pa_bnd_fill_i = model.bnd_fill_i{model.pa_bnd_i};
model.pa_bnd_fill_vec = model.bnd_fill_vec{model.pa_bnd_i};
model.pa_bnd2domain = model.bnd2domain{model.pa_bnd_i};
%
% assemble boundary condition weighting term
xnodes = model.mesh.nodes(model.pa_bnd_mesh.node_indices, 1);
ln = min(xnodes);
rn = max(xnodes);
model.beta_weight = (rn - xnodes).*(xnodes - ln)/(rn-ln)^2;
model.pa_bnd_xs = xnodes;

% assemble forcing term
force = model_opts.force_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
% weak form of the force, on nodes
if length(force) == 1
    model.b = full(sum(model.mass, 2))*force;
else
    model.b = model.mass*force(:);
end

% setup observation operator
model.obs_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.obs_locs, model.h, false);
model.n_sensors     = size(model_opts.obs_locs,1);
model.n_datasets    = 1;

% qoi
% model.qoi_mass      = model.mass_bnds{2};
model.qoi_vec       = sum(model.mass_bnds{2},1);

% constants
model.epsilon   = epsilon; % must be specified (not in the options)
model.p_rate    = p_rate;  % must be specified (not in the options)
%
model.res_tol   = model_opts.res_tol;
model.max_newton_iter       = model_opts.num_newton;
model.max_line_search_iter  = model_opts.num_linsea;

model.exp_param = model_opts.exp_param;
model.solver    = model_opts.nl_solver;
model.iter_disp = model_opts.nl_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes   = model.mesh.nodes(model.pa_bnd_mesh.node_indices, :);
n = size(nodes, 1);
D = bsxfun(@minus, nodes(:,1), nodes(:,1)').^2 + bsxfun(@minus, nodes(:,2), nodes(:,2)').^2;
D = D*inv_opts.scale;
if inv_opts.power == 2
    C = exp(-0.5*D) + 1e-10*eye(n);
else
    C = exp(-0.5*D.^(0.5*inv_opts.power));
end
prior.C = C*inv_opts.sigma^2;
prior.L = chol(prior.C, 'lower');
prior.cov_type = 'GP';
prior.type = 'Field';
prior.dof  = n;
prior.NP   = n;

if length(inv_opts.mean) == 1
    prior.mean_u = inv_opts.mean*ones(prior.dof,1);
else
    prior.mean_u = inv_opts.mean;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs.jacobian = false;
obs.like = inv_opts.like;

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
