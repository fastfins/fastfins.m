function model = build_model(model_opts, p_rate, epsilon, reg_flag)
%setup_model
%
% setup data structures used for model computation
%
% Tiangang Cui, 10/Sep/2020

if reg_flag
    [p,t] = simplemesh2d(model_opts.h, 1, model_opts.xyratio);
else
    func_mesh = @(x) drectangle0(x, 0, 1, 0, 1);
    [p, t]  = distmesh2d(func_mesh, @huniform, model_opts.h, [0, 0; 1, 1], []);
end
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
model.pa_bnd2domain = model.bnd2domain{model.pa_bnd_i};
%
% assemble boundary condition weighting term
xnodes = model.mesh.nodes(model.pa_bnd_mesh.node_indices, 1);
ln = min(xnodes);
rn = max(xnodes);
model.beta_weight = (rn - xnodes).*(xnodes - ln)/(rn-ln)^2;

% assemble forcing term
force = model_opts.force_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
% weak form of the force, on nodes
if length(force) == 1
    model.b = full(sum(model.mass, 2))*force;
else
    model.b = model.mass*force(:);
end

model.qoi_vec       = sum(model.mass_bnds{2},1);

% setup observation operator
model.obs_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.obs_locs, model.h, false);
model.n_sensors     = size(model_opts.obs_locs,1);
model.n_datasets    = 1;

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

end
