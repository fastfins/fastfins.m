function [model, obs, prior] = setup_parabolic(model_opts, inv_opts, type, output)

% Tiangang Cui, 10/Sep/2020

[p,t] = simplemesh2d(model_opts.h, 1, model_opts.xyratio);
model = setup_2nd_order(p, t, model_opts);

%
model.gmres_flag = model_opts.gmres;
model.res_tol = model_opts.res_tol;
model.explicit_ja = false;

% assemble forcing term
force = model_opts.force_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
% weak form of the force, on nodes
if length(force) == 1
    model.b = full(sum(model.mass, 2))*force;
else
    model.b = model.mass*force(:);
end

% assemble boundary condition
model.bnd_b = zeros(model.mesh.dof, 1);
model.bmass = spalloc(model.mesh.dof,model.mesh.dof,0);
%{
for bi = 1:length(model_opts.bnd_funcs)
    b_ind   = model_opts.bnd_funcs{bi}(model.mesh.nodes);
    b_force = model_opts.bc_funcs{bi}(model.mesh.nodes(:,1),model.mesh.nodes(:,2));
    if length(b_force) == 1
        b_force = ones(size(b_ind(:)))*b_force;
    else
        b_force = b_force(b_ind);
    end
    bf_weak = model.mass_bnds{bi}*sparse(b_ind(:),ones(size(b_ind(:))),b_force(:),model.mesh.dof, 1);
    switch model_opts.bc_types{bi}
        case {'flux'}
            disp('double check non-zero flux b.c.')
            model.bnd_b = model.bnd_b + bf_weak;
        case {'essential'}
            model.bmass = model.bmass + model.mass_bnds{bi}*model.penalty;
            model.bnd_b = model.bnd_b + bf_weak*model.penalty;
        case {'mixed'}
            disp('Mixed b.c. is not implemented')
    end
end
%}
model.b = model.b + model.bnd_b;


% initial condition;
model.init  = model_opts.init_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));

% apply observation operator
model.obs_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.obs_locs, model.h, false);
model.n_sensors     = size(model_opts.obs_locs,1);

model.pred_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.pred_locs, model.h, false);

%{
% apply qoi function
tol = 1E-10;
disp('flux QoI')
model.qoi_flag = false;
if ~isempty(model_opts.qoi_func)
    phi = reshape(model_opts.qoi_func(model.mesh.nodes(:,1),model.mesh.nodes(:,2)), 1, []);
    phi(abs(phi)<tol) = 0;
    model.phi = sparse(phi);
    model.qoi_flag  = true;
end
%}

model.dof = model.mesh.dof; % for consistency
% transformation
model.exp_param = model_opts.exp_param;
model.exp_thres = model_opts.exp_thres;
model.sq_param  = model_opts.sq_param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.T_obs_start   = model_opts.obs_tstart;
model.T_obs_final   = model_opts.obs_tfinal;
model.n_datasets    = model_opts.obs_ntime;
model.T_prediction  = model_opts.pred_t;
model.T_final   = max(model.T_prediction, model.T_obs_final);
model.T_lead    = model_opts.t_lead;

model.dt        = min((model.T_obs_final - model.T_obs_start)/model.n_datasets, model_opts.dt);
model.T_nstep   = ceil(model.T_final/model.dt);
model.T_obs_nstep = ceil(model.T_obs_final/model.dt);
model.obs_ind   = round(linspace(model.T_obs_start, model.T_obs_final, model.n_datasets)/model.dt);
model.pred_ind  = round(linspace(model.T_obs_final, model.T_prediction, model.n_datasets)/model.dt);

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
%
nn = 1;
switch type
    case {'vector'}
        nn = model.local_elem.num_dim;
        prior.mean_u = repmat(prior.mean_u(:), nn, 1);
    case {'tensor'}
        nn = model.local_elem.num_dim+1;
        prior.mean_u = [prior.mean_u; zeros(model.mesh.dof*(nn-1), 1)];
end
%
if nn > 1
    prior.dof = prior.dof*nn;
    switch prior.cov_type
        case {'MRF'}
            prior.M = kron(speye(nn), prior.M);
            prior.K = kron(speye(nn), prior.K);
            [prior.Lk,~,prior.pk] = chol(prior.K, 'lower', 'vector');
            [prior.Lm,~,prior.pm] = chol(prior.M, 'lower', 'vector');
        case {'GP'}
            prior.C = kron(speye(nn), prior.C);
            prior.L = kron(speye(nn), prior.L);
        case {'Conv'}
            prior.basis = kron(speye(nn), prior.basis);
            prior.basis_w = kron(speye(nn), prior.basis_w);
    end
end
prior.num_field = nn;
prior.mesh_dof = prior.mesh_dof*nn;

%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs.jacobian = false;
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
    data    = soln.d(:) + randn(size(soln.d(:)))*std;
    
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

